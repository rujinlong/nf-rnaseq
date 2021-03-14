#!/usr/bin/env nextflow
// Usage: nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "fastqc"
//        nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "clean" -resume
//        nextflow run rnaseq.nf --reads "data/*_R{1,2}.fq.gz" -profile hpc_slurm --mode "all" -resume

nextflow.enable.dsl=2

process fastqc {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/p01_fastqc"
    
    input:
    tuple val(sampleID), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: fastqc_results_ch

    when:
    params.mode == "fastqc" || params.mode == "all"

    """
    fastqc -t $task.cpus $reads
    """
}

process clean_reads {
    tag "$sampleID"
    label "big"
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

process align_to_cds {
    tag "$sampleID"
    label "big"
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


process reformat_cds_location {
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

process align_to_genome_plus_unmapped1 {
    tag "$sampleID"
    label "medium"
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

process plot_figures {
    tag "$sampleID"
    label "small"
    publishDir "$params.outdir/p05_figures"
    publishDir "$params.report/figures"
    
    input:
    tuple val(sampleID), path(depth_nonCDS), path(depth_CDS)

    output:
    path("*.png")

    when:
    params.mode == "all"

    script:
    """
    plot_figure.py -k $params.gbk_bacteria -n $depth_nonCDS -c $depth_CDS -w 2 -l 100000 -o $sampleID -f png --label
    plot_figure.py -k $params.gbk_virus -n $depth_nonCDS -c $depth_CDS -w 2 -l 20000 -o $sampleID -f png --label
    """
}

process unmapped2_reads_to_kegg {
    tag "$sampleID"
    label "medium"
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

process feature_abundance {
    label "small"
    publishDir "$params.outdir/p05_feature_counts", pattern: "*"
    publishDir "$params.report", pattern: "*"
    
    input:
    path(gff)
    path(bamfiles)

    output:
    path("*.tsv")
    path("*.gff")

    when:
    params.mode == "all"

    script:
    """
    add_non_features.py -i $gff -o with_nonCDS.gff -m $params.min_gap
    # Should be reads aligned to genome bam
    featureCounts -t CDS -g ID -a with_nonCDS.gff -o feature_counts_CDS.tsv -T $task.cpus *.bam
    # featureCounts -t gene -g ID -a with_nonCDS.gff -o feature_counts_gene.tsv -T $task.cpus *.bam
    featureCounts -t nonCDS -g ID -a with_nonCDS.gff -o feature_counts_nonCDS.tsv -T $task.cpus *.bam
    merge_feature_counts.py -c feature_counts_CDS.tsv -n feature_counts_nonCDS.tsv -o feature_counts.tsv
    """
}

process multiqc {
    label "medium"
    publishDir "$params.outdir/p98_multiqc"
    publishDir "$params.report"
    
    input:
    file(qc_rst)

    output:
    file("*")

    when:
    params.mode == "fastqc" || params.mode == "clean" || params.mode == "all"

    """
    multiqc -o . -n multiqc -s --interactive .
    """
}

workflow {
    raw_reads_ch = channel.fromFilePairs(params.reads)
    fastqc(raw_reads_ch)
    clean_reads(raw_reads_ch)
    align_to_cds(params.ref_cds, clean_reads.out.clean_reads_ch)
    reformat_cds_location(params.gbk_bacteria, params.gbk_virus, align_to_cds.out.bam_cds_ch)
    align_to_genome_plus_unmapped1(params.ref_genome, align_to_cds.out.cds_unmapped_ch.join(clean_reads.out.clean_reads_ch, by:0))
    plot_figures(align_to_genome_plus_unmapped1.out.nonCDS_depth_ch.join(reformat_cds_location.out.cds_depth_ch, by:0))
    unmapped2_reads_to_kegg(align_to_genome_plus_unmapped1.out.genome_unmapped_ch)
    feature_abundance(params.gff, align_to_genome_plus_unmapped1.out.bam_merge_ch.collect())
    
    if( params.mode == "fastqc" ) {
        multiqc(fastqc.out.fastqc_results_ch.collect())
    }
    else {
        multiqc(fastqc.out.fastqc_results_ch.concat(clean_reads.out.fastp_json_ch).concat(align_to_cds.out.bam_cds_log).concat(align_to_genome_plus_unmapped1.out.bam_genome_log).flatten().collect())
    }
}
