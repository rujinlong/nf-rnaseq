process ALIGN_TO_CDS {
    tag "$sampleID"
    publishDir "$params.outdir/p03_align_to_cds"
    
    input:
    path(ref_cds)
    tuple val(sampleID), path(reads1), path(reads2), path(singletons)

    output:
    tuple val(sampleID), path("${sampleID}_cds_unmapped.fasta.gz"), emit: cds_unmapped_ch
    path("${sampleID}_cds.log"), emit: bam_cds_log
    tuple val(sampleID), path("${sampleID}_cds.bam"), emit: bam_cds_ch

    when:
    params.mode == "all"

    script:
    """
    bowtie2-build $ref_cds reference
    
    # Mapping all clean reads to CDS
    bowtie2 --very-sensitive-local -x reference -1 $reads1 -2 $reads2 -U $singletons -p $task.cpus -S ${sampleID}_cds.sam 1> ${sampleID}_cds.log 2>&1
    samtools sort -@ $task.cpus ${sampleID}_cds.sam > ${sampleID}_cds.bam
    
    # Extract reads not mapped to CDS (will be used in process `align_to_genome_plus_unmapped1`)
    mem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
    samtools view -f 4 ${sampleID}_cds.bam | reformat.sh -Xmx\$mem in=stdin.sam out=${sampleID}_cds_unmapped.fasta.gz
    """
}