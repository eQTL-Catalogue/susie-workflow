library("SNPRelate")
library("GDSArray")
library("devtools")
library("ggplot2")
library("susieR")
library("optparse")
load_all("../eQTLUtils/")

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
  optparse::make_option(c("--outdir"), type="character", default="./finemapping_output",
                        help="Path to the output directory.", metavar = "type"),
  optparse::make_option(c("--cisdistance"), type="integer", default=1000000, 
                        help="Cis distance in bases from center of gene. [default \"%default\"]", metavar = "number"),
  optparse::make_option(c("--chunk"), type="character", default="1 1", 
                        help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Debugging
if(TRUE){
  opt = list(phenotype_list = "results/finemapping/Alasoo_2018/naive_qtl_gene_list.txt",
             cisdistance = 500000,
             gds_file = "results/finemapping/genotypes/Alasoo_2018_GRCh38.filtered.gds",
             covariates = "results/finemapping/Alasoo_2018/PCA/Alasoo_2018_ge_macrophage_naive/Alasoo_2018_ge_macrophage_naive.covariates.txt",
             expression_matrix = "results/finemapping/Alasoo_2018.gene_counts_cqn_norm.tsv",
             sample_meta = "../SampleArcheology/studies/cleaned/Alasoo_2018.tsv",
             phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz",
             chunk = "3 200",
             outdir = "./finemapping_output"
  )
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
  return(fitted)
}


#Import all files
phenotype_list = readr::read_tsv(opt$phenotype_list, col_types = "c", col_names = "phenotype_id")
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")
variant_info = importVariantInformationFromGDS(opt$gds_file)
covariates_matrix = importQtlmapCovariates(opt$covariates)

#Set parameters
cis_distance = opt$cisdistance
gds_file = opt$gds_file
study_id = sample_metadata$study[1]

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = count_matrix, 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts")

#Split phenotype list into chunks
chunk_vector = strsplit(opt$chunk, split = " ") %>% unlist() %>% as.numeric()
chunk_id = chunk_vector[1]
n_chunks = chunk_vector[2]
selected_chunk = splitIntoChunks(chunk_id, n_chunks, length(phenotype_list$phenotype_id))
selected_phenotypes = phenotype_list$phenotype_id[selected_chunk] %>% setNames(as.list(.), .)

#Split by qtl_group
qtl_groups = unique(se$qtl_group) %>% setNames(as.list(.),.)
group_se_list = purrr::map(qtl_groups, ~eQTLUtils::subsetSEByColumnValue(se, "qtl_group", .))

#Apply finemapping to all genes
results = list()
for (qtl_group in names(group_se_list)){
  message("Current QTL group is: ", qtl_group)
  results[[qtl_group]] = purrr::map(selected_phenotypes[1:2], ~finemapPhenotype(., group_se_list[[qtl_group]], variant_info, gds_file, covariates_matrix, cis_distance))
}

#Save results to disk
dir.create(opt$outdir, showWarnings = FALSE)
output_path = file.path(opt$outdir, paste0(study_id, ".", chunk_id, "_",n_chunks,".rds"))
saveRDS(results, output_path)


