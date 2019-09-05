nextflow run main.nf --cisdistance 500000 \
 --expression_matrix ~/datasets/processed/expression_matrices/normalised_studies/Alasoo_2018/Alasoo_2018.gene_counts_cqn_norm.tsv\
 --sample_meta ~/datasets/controlled_access/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
 --phenotype_meta ~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 --phenotype_list naive_qtl_gene_list.txt\
 --covariates_file ~/datasets/summary_stats/eQTLCatalogue/v0.1/pipeline_output/Alasoo_2018/PCA/Alasoo_2018_ge_macrophage_naive/Alasoo_2018_ge_macrophage_naive.covariates.txt\
 --qtl_group 'macrophage_naive'\
 --vcf ~/datasets/controlled_access/Alasoo_2018/genotypes/Alasoo_2018_GRCh38.filtered.vcf.gz\
 --n_batches 200
