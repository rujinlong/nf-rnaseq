process UNMAPPED2_READS_TO_KEGG {
    tag "$sampleID"
    publishDir "$params.outdir/p05_unmapped2_reads_to_kegg", pattern: "*.tsv"
    publishDir "$params.report", pattern: "*.tsv"
    
    input:
    tuple val(sampleID), path(genome_unmapped)

    output:
    tuple val(sampleID), path("${sampleID}_anno_KEGGvg.tsv"), emit: unmapped_vg_ch
    tuple val(sampleID), path("${sampleID}_anno_KEGGbac.tsv"), emit: unmapped_bac_ch

    when:
    params.mode == "all"

    script:
    """
    diamond blastx -p $task.cpus -d $params.kegg_vg -q $genome_unmapped -o ${sampleID}_anno_KEGGvg.tsv -f 6 --top 10 -e 0.000001
    diamond blastx -p $task.cpus -d $params.kegg_bac -q $genome_unmapped -o ${sampleID}_anno_KEGGbac.tsv -f 6 --top 10 -e 0.000001
    """
}