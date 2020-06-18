// study	qtl_group	quant_method	expression_matrix	phenotype_meta	sample_meta	vcf	phenotype_list	covariates
Channel.fromPath(params.qtl_results)
    .ifEmpty { error "Cannot find any qtl_results file in: ${params.qtl_results}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.study, row.qtl_group, row.quant_method, file(row.expression_matrix), file(row.phenotype_meta), file(row.sample_meta), file(row.vcf), file(row.phenotype_list), file(row.covariates)]}
    .set { qtl_results_ch }

process vcf_to_dosage{
    input:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates) from qtl_results_ch

    output:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates), file("${vcf.simpleName}.dose.tsv.gz"), file("${vcf.simpleName}.dose.tsv.gz.tbi") into extract_leads_ch

    script:
    if(params.vcf_genotype_field == 'DS'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nREF\\nALT\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz

        #Extract dosage and merge
        bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n" ${vcf} | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 ${vcf.simpleName}.dose.tsv.gz
        """
    } else if (params.vcf_genotype_field == 'GT'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nREF\\nALT\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz

        #Extract dosage and merge
        bcftools +dosage ${vcf} -- -t GT | tail -n+2 | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 ${vcf.simpleName}.dose.tsv.gz
        """
    } 
}

process extract_leads{
    publishDir "${params.outdir}/lead_variants/${study}/", mode: 'copy', pattern: "*.lead_variants.txt.gz"

    input:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates), file(genotype_matrix), file(genotype_matrix_index) from extract_leads_ch

    output:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file("${phenotype_list.simpleName}.lead_variants.txt.gz"), file(covariates), file(genotype_matrix), file(genotype_matrix_index) into susie_input_ch

    script:
    if(params.permuted){
        """
        zcat ${phenotype_list} | gzip > ${phenotype_list.simpleName}.lead_variants.txt.gz
        """
    } else {
        """
        zcat ${phenotype_list} | awk '{if(\$14 == 1) print \$0}' | gzip > ${phenotype_list.simpleName}.lead_variants.txt.gz
        """
    }
}

process run_susie{

    input:
    set study, qtl_group, quant_method, file(expression_matrix), file(phenotype_meta), file(sample_meta), file(vcf), file(phenotype_list), file(covariates), file(genotype_matrix), file(genotype_matrix_index) from susie_input_ch
    each batch_index from 1..params.n_batches

    output:
    set val("${study}.${qtl_group}.${quant_method}"), file("${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.txt"), file("${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.cred.txt"), file("${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}.snp.txt") into finemapping_ch

    script:
    """
    Rscript $baseDir/bin/run_susie.R --expression_matrix ${expression_matrix}\
     --phenotype_meta ${phenotype_meta}\
     --sample_meta ${sample_meta}\
     --phenotype_list ${phenotype_list}\
     --covariates ${covariates}\
     --genotype_matrix ${genotype_matrix}\
     --chunk '${batch_index} ${params.n_batches}'\
     --cisdistance ${params.cisdistance}\
     --out_prefix '${study}.${qtl_group}.${quant_method}.${batch_index}_${params.n_batches}'\
     --qtl_group ${qtl_group}\
     --eqtlutils ${params.eqtlutils}\
     --permuted ${params.permuted}
    """
}

process merge_susie{
    publishDir "${params.outdir}/susie/", mode: 'copy', pattern: "*.txt.gz"

    input:
    set study_qtl_group_quant, in_cs_variant_batch_names, credible_set_batch_names, variant_batch_names from finemapping_ch.groupTuple()

    output:
    file("${study_qtl_group_quant}.txt.gz"), file("${study_qtl_group_quant}.cred.txt.gz"), file("${study_qtl_group_quant}.snp.txt.gz") into susie_merged_ch

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${in_cs_variant_batch_names.join(' ')} | bgzip -c > ${study_qtl_group_quant}.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${credible_set_batch_names.join(' ')} | bgzip -c > ${study_qtl_group_quant}.cred.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${variant_batch_names.join(' ')} | bgzip -c > ${study_qtl_group_quant}.snp.txt.gz
    """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
}

