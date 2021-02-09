nextflow.enable.dsl=2

// Input channel for most finemapping inputs
Channel.fromPath(params.qtl_results)
    .ifEmpty { error "Cannot find any qtl_results file in: ${params.qtl_results}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.count_matrix), file(row.pheno_meta), file(row.sample_meta), file(row.phenotype_list), file(row.covariates)]}
    .set { qtl_results_ch }

// Separate input channel for the VCF file
Channel.fromPath(params.qtl_results)
    .ifEmpty { error "Cannot find any qtl_results file in: ${params.qtl_results}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.vcf) ]}
    .set { vcf_ch }

batch_ch = Channel.of(1..params.n_batches)

include { vcf_to_dosage } from './modules/vcf_to_dosage'
include { run_susie; merge_susie; sort_susie } from './modules/susie'

workflow {
    vcf_to_dosage(vcf_ch)
    run_susie(qtl_results_ch.join(vcf_to_dosage.out), batch_ch)
    merge_susie( run_susie.out.groupTuple(size: params.n_batches) )
    sort_susie( merge_susie.out )
}