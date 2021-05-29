process CLEAN_READS {
    tag "$sampleID"
    publishDir "$params.outdir/p02_reads_clean", pattern: "*.gz"
    
    input:
    tuple val(sampleID), path(reads)

    output:
    path("*.json"), emit: fastp_json_ch
    tuple val(sampleID), path("${sampleID}_clean_R1.fq.gz"), path("${sampleID}_clean_R2.fq.gz"), path("${sampleID}_singletons.fq.gz"), emit: clean_reads_ch

    when:
    params.mode == "clean" || params.mode == "all"

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sampleID}_clean_R1.fq.gz -O ${sampleID}_clean_R2.fq.gz --unpaired1 ${sampleID}_singletons.fq.gz --unpaired2 ${sampleID}_singletons.fq.gz --failed_out ${sampleID}_fail.fq.gz -f $params.leftcutR1 -t $params.rightcutR1 -F $params.leftcutR2 -T $params.rightcutR2 --detect_adapter_for_pe -p -w $task.cpus -n 1 -l 30 -5 -W 4 -M 20 -r -c -g -x -j ${sampleID}_fastp.json -h ${sampleID}_fastp.html
    """
}