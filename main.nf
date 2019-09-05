Channel
    .fromPath(params.vcf)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.vcf}" } 
    .set { genotype_vcf }

Channel
    .fromPath(params.expression_matrix)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.expression_matrix}" } 
    .set { expression_matrix_ch }

Channel
    .fromPath(params.phenotype_meta)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.phenotype_meta}" } 
    .set { phenotype_meta_ch }

Channel
    .fromPath(params.sample_meta)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.sample_meta}" } 
    .set { sample_meta_ch }

Channel
    .fromPath(params.phenotype_list)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.phenotype_list}" } 
    .set { phenotype_list_ch }

Channel
    .fromPath(params.covariates_file)
    .ifEmpty { exit 1, "VCF genotype file not found: ${params.covariates}" } 
    .set { covariates_file_ch }

process vcf_to_gds{

    input:
    file input_vcf from genotype_vcf

    output:
    file "study_genotypes.gds" into genotype_gds_ch

    script:
    """
    Rscript $baseDir/bin/vcf_to_gds.R --vcf ${input_vcf} --gds study_genotypes.gds
    """
}

process run_susie{

    input:
    file gds_file from genotype_gds_ch
    file expression_matrix from expression_matrix_ch
    file phenotype_meta from phenotype_meta_ch
    file sample_meta from sample_meta_ch
    file phenotype_list from phenotype_list_ch
    file covariates_file from covariates_file_ch
    each batch_index from 1..params.n_batches

    output:
    file("finemapping_results.${batch_index}_${params.n_batches}.txt") into finemapping_ch

    script:
    """
    Rscript $baseDir/bin/run_susie.R --expression_matrix ${expression_matrix}\
     --phenotype_meta ${phenotype_meta}\
     --sample_meta ${sample_meta}\
     --phenotype_list ${phenotype_list}\
     --covariates ${covariates_file}\
     --gds_file ${gds_file}\
     --chunk '${batch_index} ${params.n_batches}'\
     --cisdistance ${params.cisdistance}\
     --outfile 'finemapping_results.${batch_index}_${params.n_batches}.txt'\
     --qtl_group ${params.qtl_group}\
     --eqtlutils ${params.eqtlutils}
    """
}

