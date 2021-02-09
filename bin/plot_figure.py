#!/usr/bin/env python

import click
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

def remove_label(x):
    x.label = None
    return x
        
def plot_figure(record, depth, window_size, line_len, fout, label=False):
    """
    input:
    ------
    rec: 
    ncdepth: non-CDS depth
    gdepth: genome depth
    fout: saved figure file name
    """
    genome_len = int(depth.location.iloc[-1])
    ymax = max(depth.nonCDS_log2.max(), depth.CDS_log2.max())*2
    nfig = int(genome_len/line_len)+1
    fig = plt.figure(figsize=(100, 7*nfig))
    outer = gridspec.GridSpec(nfig, 1, wspace=0, hspace=0)
    for i in range(0, nfig):
        start = i*line_len
        end = start + line_len
        if end >= genome_len:
            end = genome_len
        rec = record.crop([start, end])
        
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, 
                                                 subplot_spec=outer[i], 
                                                 wspace=0, 
                                                 hspace=0,
                                                 height_ratios=[5,2])
    
        ## Genome annotation plot
        if not label:
            rec.features = [remove_label(x) for x in rec.features]
        ax1 = plt.Subplot(fig, inner[0])
        rec.plot(ax=ax1, with_ruler=True, strand_in_label_threshold=2)
        fig.add_subplot(ax1)

        ## Mapping to non-CDS region
        ax2 = plt.Subplot(fig, inner[1])
        depth.loc[start:end, ["CDS_log2", "nonCDS_log2"]].rolling(window=window_size, center=True).mean().fillna(0).plot.area(ax=ax2, alpha=0.3, color=["green", "red"])
        ax2.set_ylabel("#Reads (log2)")
        ax2.set_ylim(0,ymax)
        fig.add_subplot(ax2)

    ## Save figure
    fig.savefig(fout)


@click.command()
@click.option("--fname_gbk", '-k', help="Reference genbank file")
@click.option("--fname_noncds_depth", '-n', help="Depth file of non CDS mapping")
@click.option("--fname_cds_depth", '-c', help="Depth file of cds mapping")
@click.option("--window_size", '-w', type=int, help="Window size")
@click.option("--line_len", '-l', type=int, help="Length of each line")
@click.option("--figfmt", '-f', help="Figure format [png, pdf]")
@click.option("--fout", '-o', help="output file name prefix")
@click.option('--label/--no-label', default=False, help="Plot label or not")
def main(fname_gbk, fname_noncds_depth, fname_cds_depth, window_size, line_len, figfmt, fout, label):
    gbk = SeqIO.read(fname_gbk, "genbank")
    gbk.id = gbk.name
    ndepth = pd.read_csv(fname_noncds_depth, sep='\t', names=['ctgid', 'location', 'nonCDS'], skiprows=1)
    ndepth = ndepth[ndepth.ctgid==gbk.id].reset_index(drop=True)
    cdepth = pd.read_csv(fname_cds_depth, sep='\t', names=['ctgid', 'location', 'CDS'], skiprows=1)
    cdepth = cdepth[cdepth.ctgid==gbk.id].reset_index(drop=True)
   
    df = cdepth.merge(ndepth, on=['ctgid', 'location'], how="outer")
    df.fillna(0, inplace=True)
    df['nonCDS_log2'] = np.log2(df['nonCDS']+1)
    df['CDS_log2'] = np.log2(df['CDS']+1)

    color_map = {gbk.id: "blue", "CDS": "yellow"}
    rec = BiopythonTranslator(features_filters=(lambda f: f.type in ["CDS", gbk.id],),
                              features_properties=lambda f: {"color": color_map.get(f.type, "white")}).translate_record(gbk)
    ## Add genome feature
    # rec.features+=rec2.features
    if len(df) > 0:
        plot_figure(rec, df, window_size, line_len, '{}_{}.{}'.format(fout, gbk.id, figfmt), label)


if __name__ == '__main__':
    main()
