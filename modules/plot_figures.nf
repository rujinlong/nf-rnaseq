process PLOT_FIGURES {
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