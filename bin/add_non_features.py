#!/usr/bin/env python

import click
import pandas as pd

def add_non_features(df, clms, min_gap=50):
    df.sort_values(by='nstart', inplace=True)
    nprev = 1
    non_features = []
    for i in range(len(df)):
        nstart = df.loc[i, 'nstart']
        if nstart - nprev > min_gap:
            new_start = nprev
            new_end = nstart
            seqname = df.loc[i, 'seqname']
            non_features.append([seqname, 
                               df.loc[i, 'source'],
                               'nonCDS',
                               nprev,
                               nstart,
                               '.',
                               df.loc[i, 'strand'],
                               df.loc[i, 'frame'],
                               "ID={}:{}-{}_nonCDS".format(seqname, str(nprev), str(nstart))])
        nprev = df.loc[i, 'nend']

    df_non_features = pd.DataFrame(non_features, columns=clms)
    tbl = pd.concat([df, df_non_features], axis=0)
    tbl.sort_values(by=['seqname', 'nstart'], inplace=True)
    tbl.reset_index(drop=True, inplace=True)
    return tbl


@click.command()
@click.option("--fin", '-i', help="input GFF file name")
@click.option("--fout", '-o', help="output GFF file name")
@click.option("--min_gap", '-m', default=50, type=int, help="Minimum gap between two features.")
def main(fin, fout, min_gap):
    """
    Usage:
        add_non_features.py -i hostphage.gff -o hostphage_with_nonfeature.gff -m 50
    """
    clms = ['seqname', 'source', 'feature', 'nstart', 'nend', 'score', 'strand', 'frame', 'attr']
    df = pd.read_csv(fin, sep='\t', comment='#', names=clms)
    # Use CDS, instead of gene, because no gene annotation for phage
    df = df[df.feature=='CDS'].copy()
    df.reset_index(drop=True, inplace=True)
    print(df.shape)
    df2 = add_non_features(df, clms, min_gap)
    df2.to_csv(fout, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()

