#!/usr/bin/env python

import click
import pandas as pd


@click.command()
@click.option("--fcds", '-c', help="feature_counts_CDS.tsv")
@click.option("--fnoncds", '-n', help="feature_counts_nonCDS.tsv")
@click.option("--fout", '-o', help="feature_counts_all.tsv")
def main(fcds, fnoncds, fout):
    cds = pd.read_csv(fcds, sep='\t', comment='#')
    noncds = pd.read_csv(fnoncds, sep='\t', comment='#')
    
    clms1 = cds.columns[:6].to_list()
    samples = cds.columns[6:].to_list()
    samples.sort()
    clms = clms1 + samples
    
    df = pd.concat([cds[clms], noncds[clms]], axis=0)
    df.sort_values(by=['Chr', 'Start'], inplace=True)
    df.columns = [x.replace('_genome.bam', '') for x in df.columns]
    df.to_csv(fout, sep='\t', index=False)


if __name__ == '__main__':
    main()