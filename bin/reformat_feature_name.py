#!/usr/bin/env python

import click
import pandas as pd
import re

def create_feature_name(x):
    if x['Geneid'].startswith('vB'):
        new_name = 'v_{}-{}'.format(x['Start'], x['End'])
    else:
        new_name = '{}_{}-{}'.format(re.sub(r'^cds-', '', x['Geneid']), x['Start'], x['End'])
    return new_name


def reformat_table(fin, fout):
    # fin = 'rnaseqd0175/feature_counts_cds.tsv'
    df = pd.read_csv(fin, sep='\t', comment="#")
    df.columns = [re.sub(r'_genome.bam$', '', x) for x in df.columns.tolist()]
    sample_ids = df.columns.tolist()[6:]
    meta_clm = df.columns.tolist()[:6]
    
    df['feature_name'] = df.apply(lambda x:create_feature_name(x), axis=1)
    df[['feature_name']+meta_clm+sample_ids].to_csv(fout, sep='\t', index=False)
    


@click.command()
@click.option("--fin", '-i', help="input file name")
@click.option("--fout", '-o', default='out.txt', help="output file name")
def main(fin, fout):
    reformat_table(fin, fout)


if __name__ == '__main__':
    main()
