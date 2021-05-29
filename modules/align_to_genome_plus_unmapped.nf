process ALIGN_TO_GENOME_PLUS_UNMAPPED1 {
    tag "$sampleID"
    publishDir "$params.outdir/p04_align_to_genome"
    publishDir "$params.report/depth",  pattern: "${sampleID}_nonCDS_genome.depth"
    
    input:
    path(ref_genome)
    tuple val(sampleID), path(cds_unmapped), path(reads1), path(reads2), path(singletons)

    output:
    path("*.depth")
    path("*_genome.log"), emit: bam_genome_log
    path("${sampleID}_genome.bam"), emit: bam_merge_ch
    tuple val(sampleID), path("${sampleID}_genome_unmapped.fasta.gz"), emit: genome_unmapped_ch
    tuple val(sampleID), path("${sampleID}_nonCDS_genome.depth"), emit: nonCDS_depth_ch

    when:
    params.mode == "all"

    script:
    """
    bowtie2-build $ref_genome reference
    
    # Mapping all clean reads to genome
    bowtie2 --very-sensitive-local -x reference -1 $reads1 -2 $reads2 -U $singletons -p $task.cpus -S ${sampleID}_genome.sam 1> ${sampleID}_genome.log 2>&1
    samtools sort -@ $task.cpus ${sampleID}_genome.sam > ${sampleID}_genome.bam
    samtools depth ${sampleID}_genome.bam > tmp_genome.depth
    calc_depth.py -i tmp_genome.depth -o ${sampleID}_genome.depth
    
    # Extract unmapped reads (will be used in process `unmapped2_reads_to_kegg` )
    mem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
    samtools view -f 4 ${sampleID}_genome.bam | reformat.sh -Xmx\$mem in=stdin.sam out=${sampleID}_genome_unmapped.fasta.gz

    # Mapping non-CDS reads to genome
    bowtie2 --very-sensitive-local -x reference -f $cds_unmapped -p $task.cpus -S ${sampleID}_nonCDS_genome.sam 1> ${sampleID}_nonCDS_genome.log 2>&1
    samtools sort -@ $task.cpus ${sampleID}_nonCDS_genome.sam > ${sampleID}_nonCDS_genome.bam
    samtools depth ${sampleID}_nonCDS_genome.bam > tmp_nonCDS.depth
    calc_depth.py -i tmp_nonCDS.depth -o ${sampleID}_nonCDS_genome.depth
    rm tmp_*.depth
    """
}