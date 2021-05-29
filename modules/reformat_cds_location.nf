process REFORMAT_CDS_LOCATION {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/p04_reformat_cds_location"
    publishDir "$params.report/depth",  pattern: "${sampleID}_cds.depth"
    
    input:
    path(gbk_bacteria)
    path(gbk_virus)
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("${sampleID}_cds.depth"), emit: cds_depth_ch

    when:
    params.mode == "all"

    script:
    """
    samtools depth $bam > tmp.depth
    reformat_cds_location.py -g $gbk_bacteria -d tmp.depth -t bacteria -i $params.bacteria_id -o ${sampleID}_cds_bacteria.depth
    reformat_cds_location.py -g $gbk_virus -d tmp.depth -t virus -i $params.virus_id -o ${sampleID}_cds_virus.depth
    cat ${sampleID}_cds_*.depth >> ${sampleID}_merge.depth
    calc_depth.py -i ${sampleID}_merge.depth -o ${sampleID}_cds.depth
    rm ${sampleID}_cds_*.depth ${sampleID}_merge.depth tmp.depth
    """
}