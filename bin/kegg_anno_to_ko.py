#!/usr/bin/env python

import click
import pandas as pd

@click.command()
@click.option("--fname_kegg_anno", '-a', help="Results from diamond search against KEGG protein database")
@click.option("--fname_gene2ko", '-m', help="gene2ko_nr.map")
@click.option("--out", '-o', help="protein_to_ko.tsv")
@click.option("--bit_score", '-b', default=50, type=float, help="bit score")
@click.option("--evalue", '-e', default=0.001, type=float, help="evalue")
def main(fname_kegg_anno, fname_gene2ko, out, bit_score=50, evalue=0.001):
    clms = ['prot_id', 'kegg_gene_id', 'pct_identity', 'ali_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
    df_anno = pd.read_csv(fname_kegg_anno, sep='\t', names=clms)
    df_anno_sel = df_anno[(df_anno.bit_score > bit_score) & (df_anno.evalue < evalue)]
    df_left = df_anno_sel[['prot_id', 'kegg_gene_id']].copy()
    
    df_g2ko = pd.read_csv(fname_gene2ko, sep='\s+', names=['kegg_gene_id', 'koid'])
    df_final = df_left.merge(df_g2ko, on='kegg_gene_id', how='inner')
    df_final.to_csv(out, sep='\t', index=False)
    
if __name__ == '__main__':
    main()
