#!/usr/bin/env nextflow
// Usage: nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "fastqc"
//        nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "clean" -resume
//        nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "all" -resume

nextflow.enable.dsl=2


log.info """\
NF - RNASEQ PIPELINE
=========================
result: ${params.outdir}
report: ${params.report}
"""


include { FASTQC } from './modules/fastqc'
include { CLEAN_READS } from './modules/clean_reads'
include { ALIGN_TO_CDS } from './modules/align_to_cds'
include { REFORMAT_CDS_LOCATION } from './modules/reformat_cds_location'
include { ALIGN_TO_GENOME_PLUS_UNMAPPED1 } from './modules/align_to_genome_plus_unmapped'
include { PLOT_FIGURES } from './modules/plot_figures'
include { UNMAPPED2_READS_TO_KEGG } from './modules/unmapped2_reads_to_kegg'
include { FEATURE_ABUNDANCE } from './modules/feature_abundance'
include { PATHWAY_ANALYSIS } from './modules/pathway_analysis'
include { CLEAN_OUTPUT } from './modules/clean_output'
include { MULTIQC } from './modules/multiqc'


workflow {
    raw_reads_ch = channel.fromFilePairs(params.reads)
    FASTQC(raw_reads_ch)
    CLEAN_READS(raw_reads_ch)
    ALIGN_TO_CDS(params.ref_cds, CLEAN_READS.out.clean_reads_ch)
    REFORMAT_CDS_LOCATION(params.gbk_bacteria, params.gbk_virus, ALIGN_TO_CDS.out.bam_cds_ch)
    ALIGN_TO_GENOME_PLUS_UNMAPPED1(params.ref_genome, ALIGN_TO_CDS.out.cds_unmapped_ch.join(CLEAN_READS.out.clean_reads_ch, by:0))
    PLOT_FIGURES(ALIGN_TO_GENOME_PLUS_UNMAPPED1.out.nonCDS_depth_ch.join(REFORMAT_CDS_LOCATION.out.cds_depth_ch, by:0))
    UNMAPPED2_READS_TO_KEGG(ALIGN_TO_GENOME_PLUS_UNMAPPED1.out.genome_unmapped_ch)
    FEATURE_ABUNDANCE(params.gff, ALIGN_TO_GENOME_PLUS_UNMAPPED1.out.bam_merge_ch.collect())
    PATHWAY_ANALYSIS(params.viral_cds)
    CLEAN_OUTPUT(FEATURE_ABUNDANCE.out.feature_counts, PATHWAY_ANALYSIS.out.cds2ko)
    
    if( params.mode == "fastqc" ) {
        MULTIQC(FASTQC.out.fastqc_results_ch.collect())
    }
    else {
        MULTIQC(FASTQC.out.fastqc_results_ch.concat(CLEAN_READS.out.fastp_json_ch).concat(ALIGN_TO_CDS.out.bam_cds_log).concat(ALIGN_TO_GENOME_PLUS_UNMAPPED1.out.bam_genome_log).flatten().collect())
    }
}
