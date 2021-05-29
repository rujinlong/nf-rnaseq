process PATHWAY_ANALYSIS {
    publishDir "$params.outdir/p06_pathway", pattern: "*"
    publishDir "$params.report", pattern: "*"
    
    input:
    path(viral_cds)

    output:
    path("*.tsv")
    path("cds2ko.tsv"), emit: cds2ko

    when:
    params.mode == "all"

    script:
    """
    diamond blastx -p $task.cpus -d $params.kegg_db -q $viral_cds -o anno_virus_KEGGnvg.diamond -f 6 --top 1 -e 0.000001
    diamond blastx -p $task.cpus -d $params.kegg_vg -q $viral_cds -o anno_virus_KEGGvg.diamond -f 6 --top 1 -e 0.000001
    cat anno_virus_KEGG* > anno_virus_mergeKEGG.diamond
    kegg_anno_to_ko.py -a anno_virus_mergeKEGG.diamond -m $params.kegg_gene2ko -o tmp.tsv -b 50 -e 0.000001
    kegg_anno_to_ko.py -a anno_virus_mergeKEGG.diamond -m $params.keggvg_gene2ko -o tmp2.tsv -b 50 -e 0.000001
    cat tmp.tsv tmp2.tsv | grep -v '^prot_id' | cut -f1,3 | sort -u > virgene2ko.tsv
    cat $params.bacgene2ko virgene2ko.tsv > cds2ko.tsv
    """
}