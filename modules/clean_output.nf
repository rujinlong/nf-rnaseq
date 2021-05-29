process CLEAN_OUTPUT {
    label "small"
    publishDir "$params.outdir/p07_clean_output", pattern: "*"
    publishDir "$params.report", pattern: "*"
    
    input:
    path(feature_counts)
    path(cds2ko)

    output:
    path("*.tsv")

    when:
    params.mode == "all"

    script:
    """
    kegg_pathway.py -b $params.gff_bacteria -v $params.gbk_virus -c $feature_counts -k $cds2ko -p $params.kegg_pathway -o table_feature_counts.tsv -f table_annotation_counts.tsv
    """
}