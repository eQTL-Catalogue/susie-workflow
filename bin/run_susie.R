suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("GDSArray"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("susieR"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--phenotype_meta"), type="character", default=NULL,
                        help="Phenotype metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--sample_meta"), type="character", default=NULL,
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("--expression_matrix"), type="character", default=NULL,
                        help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("--phenotype_list"), type="character", default=NULL,
                        help="Path to the phenotype list file.", metavar = "type"),
  optparse::make_option(c("--gds_file"), type="character", default=NULL,
                        help="Raw genotypes in gds format.", metavar = "type"),
  optparse::make_option(c("--covariates"), type="character", default=NULL,
                        help="Path to covariates file in QTLtools format.", metavar = "type"),
  optparse::make_option(c("--outtxt"), type="character", default="./finemapping_output.txt",
                        help="Name of the output .txt file.", metavar = "type"),
  optparse::make_option(c("--qtl_group"), type="character", default=NULL,
                        help="Value of the current qtl_group.", metavar = "type"),
  optparse::make_option(c("--cisdistance"), type="integer", default=1000000, 
                        help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("--chunk"), type="character", default="1 1", 
                        help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  optparse::make_option(c("--eqtlutils"), type="character", default=NULL,
              help="Optional path to the eQTLUtils R package location. If not specified then eQTLUtils is assumed to be installed in the container. [default \"%default\"]", metavar = "type"),
  optparse::make_option(c("--permuted"), type="character", default="true",
                        help="If 'false', lead variants were extracted from nominal p-value files. [default \"%default\"]", metavar = "type")
  
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Debugging
if(FALSE){
  opt = list(phenotype_list = "testdata/monocyte_LPS2.permuted.txt.gz",
             cisdistance = 200000,
             gds_file = "testdata/Fairfax_2014_GRCh38.filtered.renamed.gds",
             covariates = "testdata/monocyte_LPS2.covariates.txt",
             expression_matrix = "testdata/Fairfax_2014.HumanHT-12_V4_norm_exprs.tsv.gz",
             sample_meta = "testdata/Fairfax_2014.tsv",
             phenotype_meta = "testdata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz",
             chunk = "3 200",
             outtxt = "./finemapping_output.txt",
             eqtlutils = "../eQTLUtils/",
             qtl_group = "monocyte_LPS2",
             permuted = "true"
  )
}

#Load eQTLUtils
if(opt$eqtlutils == "null"){
  opt$eqtlutils = NULL
}
if (!is.null(opt$eqtlutils)){
  devtools::load_all(opt$eqtlutils)
}

#Print all options
print(opt)

#Define helper functions
importQtlmapCovariates <- function(covariates_path){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$SampleID
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  return(pc_matrix)
}

splitIntoBatches <- function(n, batch_size){
  n_batches = ceiling(n/batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total){
  chunk_size = floor(n_total/(n_chunks))
  batches = splitIntoBatches(n_total,chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

finemapPhenotype <- function(phenotype_id, se, variant_info, gds_file, covariates, cis_distance){
  message("Processing phenotype: ", phenotype_id)
  
  #Extract phenotype from SE
  gene_vector = eQTLUtils::extractPhentypeFromSE(phenotype_id, se, "counts") %>%
    dplyr::mutate(phenotype_value_std = (phenotype_value - mean(phenotype_value))/sd(phenotype_value))
  selected_phenotype = phenotype_id
  gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == selected_phenotype)
  
  #Rearrange samples in the covariates matrix
  covariates_matrix = covariates[gene_vector$genotype_id,]
  
  #Import genotype matrix
  genotype_matrix = eQTLUtils::extractGenotypeMatrixFromGDS(
    chr = gene_meta$chromosome, 
    start = gene_meta$phenotype_pos - cis_distance, 
    end = gene_meta$phenotype_pos + cis_distance, 
    variant_information = variant_info, 
    gdsfile = gds_file)
  
  #Residualise gene expression
  model_fit = stats::lm.fit(covariates_matrix, gene_vector$phenotype_value_std)
  residuals = dplyr::mutate(gene_vector, phenotype_residual = model_fit$residuals) %>%
    dplyr::mutate(phenotype_residual_std = (phenotype_residual - mean(phenotype_residual))/sd(phenotype_residual))
  
  #Fit finemapping model
  expression_vector = residuals$phenotype_residual_std
  names(expression_vector) = residuals$genotype_id
  gt_matrix = genotype_matrix[,names(expression_vector)]
  gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))
  
  fitted <- susieR::susie(gt_std, expression_vector,
                          L = 10,
                          estimate_residual_variance = TRUE, 
                          estimate_prior_variance = FALSE,
                          scaled_prior_variance = 0.1,
                          verbose = TRUE,
                          compute_univariate_zscore = TRUE)
  fitted$variant_id = rownames(gt_matrix)
  return(fitted)
}

extractCredibleSets <- function(susie_object){
  credible_sets = susie_object$sets$cs
  cs_list = list()
  for (index in seq_along(credible_sets)){
    cs_variants = credible_sets[[index]]
    cs_list[[index]] = dplyr::data_frame(cs_id = paste0("L", index),
                                         pip = susie_object$pip[cs_variants],
                                         variant_id = susie_object$variant_id[cs_variants],
                                         z = susie_object$z[cs_variants],
                                         converged = susie_object$converged)
  }
  df = purrr::map_df(cs_list, identity)

  #Extract purity values for all sets
  purity_res = susie_object$sets$purity
  set_ids = rownames(purity_res)
  purity_df = dplyr::as_tibble(purity_res) %>%
    dplyr::mutate(cs_id = rownames(purity_res))

  if(nrow(df) > 0 & nrow(purity_df) > 0){
    res_df = dplyr::left_join(df, purity_df, by = "cs_id")
  } else{
    res_df = df
  }
  return(res_df)
}


#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")
variant_info = eQTLUtils::importVariantInformationFromGDS(opt$gds_file)
covariates_matrix = importQtlmapCovariates(opt$covariates)

#Import list of phenotypes for finemapping
if (opt$permuted == "true"){
  phenotype_table = eQTLUtils::importQTLtoolsTable(opt$phenotype_list)
  filtered_list = dplyr::filter(phenotype_table, p_fdr < 0.01)
  phenotype_list = dplyr::semi_join(phenotype_meta, filtered_list, by = "group_id")
  message("Number of phenotypes included for analysis: ", nrow(phenotype_list))
} else {
  nominal_leads = read.table(opt$phenotype_list, stringsAsFactors = FALSE) %>%
    dplyr::transmute(phenotype_id = V1, n_variants = V6, p_nominal= V12) %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(phenotype_id) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, n = n_variants)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>%
    dplyr::filter(p_fdr < 0.1) #More lenient threshold for bonferroni correction
  phenotype_list = dplyr::semi_join(phenotype_meta, nominal_leads, by = "phenotype_id")
  message("Number of phenotypes included for analysis: ", nrow(phenotype_list))
}

#Set parameters
cis_distance = opt$cisdistance
gds_file = opt$gds_file
study_id = sample_metadata$study[1]

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts",
                                                         reformat = FALSE)

#Split phenotype list into chunks
chunk_vector = strsplit(opt$chunk, split = " ") %>% unlist() %>% as.numeric()
chunk_id = chunk_vector[1]
n_chunks = chunk_vector[2]
selected_chunk = splitIntoChunks(chunk_id, n_chunks, length(phenotype_list$phenotype_id))
selected_phenotypes = phenotype_list$phenotype_id[selected_chunk] %>% setNames(as.list(.), .)

#Check that the qtl_group is valid and subset
assertthat::assert_that(opt$qtl_group %in% unique(se$qtl_group))
selected_qtl_group = eQTLUtils::subsetSEByColumnValue(se, "qtl_group", opt$qtl_group)

#Apply finemapping to all genes
results = purrr::map(selected_phenotypes, ~finemapPhenotype(., selected_qtl_group, 
                                                            variant_info, gds_file, covariates_matrix, cis_distance))

#Extract credible sets from finemapping results
cs_df = purrr::map_df(results, extractCredibleSets, .id = "phenotype_id") 
if(nrow(cs_df) > 0){
  cs_df = dplyr::group_by(phenotype_id, cs_id) %>%
    dplyr::mutate(cs_size = n()) %>%
    dplyr::ungroup() %>%
    tidyr::separate(variant_id, c("chr", "pos", "ref", "alt"),sep = "_", remove = FALSE) %>%
    dplyr::mutate(chr = stringr::str_remove_all(chr, "chr")) %>%
    dplyr::transmute(phenotype_id, variant_id, cs_id, chr, pos, ref, alt, pip, z, cs_min_r2 = min.abs.corr, cs_avg_r2 = mean.abs.corr, cs_size)
}
write.table(cs_df, opt$outtxt, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


