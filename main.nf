// study	qtl_group	quant_method	expression_matrix	phenotype_meta	sample_meta	vcf	phenotype_list	covariates
Channel.fromPath(params.qtl_results)
    .ifEmpty { error "Cannot find any qtl_results file in: ${params.qtl_results}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.study, row.qtl_group, row.quant_method, file(row.expression_matrix), file(row.phenotype_meta), file(row.sample_meta), file(row.vcf), file(row.phenotype_list), file(row.covariates)]}
    .set { qtl_results_ch }

process vcf_to_gds{

    input:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates) from qtl_results_ch

    output:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates), file("${vcf.simpleName}.gds") into susie_input_ch

    script:
    """
    Rscript $baseDir/bin/vcf_to_gds.R --vcf ${vcf} --gds ${vcf.simpleName}.gds
    """
}

process run_susie{
    publishDir "${params.outdir}/susie/${study}", mode: 'copy'


    input:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates), file(gds) from susie_input_ch
    each batch_index from 1..params.n_batches

    output:
    set file("${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.rds"), file("${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.txt") into finemapping_ch

    script:
    """
    Rscript $baseDir/bin/run_susie.R --expression_matrix ${expression_matrix}\
     --phenotype_meta ${phenotype_meta}\
     --sample_meta ${sample_meta}\
     --phenotype_list ${phenotype_list}\
     --covariates ${covariates}\
     --gds_file ${gds}\
     --chunk '${batch_index} ${params.n_batches}'\
     --cisdistance ${params.cisdistance}\
     --outrds '${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.rds'\
     --outtxt '${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.txt'\
     --qtl_group ${qtl_group}\
     --eqtlutils ${params.eqtlutils}
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}

