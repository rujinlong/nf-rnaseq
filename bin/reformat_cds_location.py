#!/usr/bin/env python

import click
import pandas as pd
from Bio import SeqIO


@click.command()
@click.option("--fgbk", '-g', help="Genome in genbank format")
@click.option("--fdepth", '-d', help="CDS depth file")
@click.option("--taxa", '-t', help="[bacteria, virus]")
@click.option("--genome_id", '-i', help="[bacteria, virus]")
@click.option("--fout", '-o', help="output file name")
def main(fgbk, fdepth, taxa, genome_id, fout):
    """
    taxa: only works for vB_AB_HMGU1
    """
    gbk = SeqIO.read(fgbk, 'genbank')
    depth = pd.read_csv(fdepth, sep='\t', names=['cdsid', 'location', 'depth'])
    depth['ctgid'] = depth.apply(lambda x:x['cdsid'].split(':')[0], axis=1)
    depth = depth[depth.ctgid==genome_id]
    if len(depth) > 0:
        if taxa == "bacteria":
            features = [[f.qualifiers['protein_id'][0], f.location.start.position] for f in gbk.features if f.type=='CDS']
            depth['feature_id'] = depth.apply(lambda x: x['cdsid'].split(':')[1].split('_')[0], axis=1)
        elif taxa == "virus":
            features = [["{}-{}".format(f.location.start.position, f.location.end.position), f.location.start.position] for f in gbk.features]
            depth['feature_id'] = depth.apply(lambda x: x['cdsid'].split(':')[1], axis=1)
        features_start = pd.DataFrame(features, columns=['feature_id', 'start_location'])
    
        df = depth.merge(features_start, on='feature_id', how='left')
        df['location'] = df['location'] + df['start_location']-1
        df[['ctgid', 'location', 'depth']].to_csv(fout, index=False, header=False, sep='\t')


if __name__ == '__main__':
    main()