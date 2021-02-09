#!/usr/bin/env python

import click
import pandas as pd
import numpy as np

def depth_of_contig_locations(incomplete_depth_locations):
    contig_len = incomplete_depth_locations[-1][0]
    loc_depth = np.array([0] * contig_len)
    for i in range(len(incomplete_depth_locations)):
        idx = incomplete_depth_locations[i][0] - 1
        loc_depth[idx] = incomplete_depth_locations[i][1]

    return loc_depth
    
    
def depth_of_all_contigs(depth_report):
    df = pd.read_csv(depth_report, sep='\t', names=['ctgid', 'location', 'depth'])
    contigs = df.ctgid.unique()
    rst = []
    for ctg in contigs:
        incomplete_depth_locations = df[df.ctgid==ctg][['location', 'depth']].values
        if len(incomplete_depth_locations) > 0:
            complete_depth_locations = depth_of_contig_locations(incomplete_depth_locations)
            tmp = pd.DataFrame(complete_depth_locations, columns=['depth'])
            tmp['ctgid'] = ctg
            tmp['location'] = tmp.index+1
            rst.append(tmp[['ctgid', 'location', 'depth']])
    tbl = pd.concat(rst, axis=0)
    return tbl


@click.command()
@click.option("--fin", '-i', help="depth file")
@click.option("--fout", '-o', help="output file name")
def main(fin, fout):
    df = depth_of_all_contigs(fin)
    df.to_csv(fout, index=False, sep='\t')


if __name__ == '__main__':
    main()