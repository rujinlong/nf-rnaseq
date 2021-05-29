process FEATURE_ABUNDANCE {
    label "small"
    publishDir "$params.outdir/p05_feature_counts", pattern: "*"
    publishDir "$params.report", pattern: "*"
    
    input:
    path(gff)
    path(bamfiles)

    output:
    path("*.tsv")
    path("*.gff")
    path("feature_counts.tsv"), emit: feature_counts

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