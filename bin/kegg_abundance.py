
import pandas as pd
from argparse import ArgumentParser


def kegg_output(tbl, fid):
    anno_col = 'KEGG_Annotation'
    df = pd.read_csv(tbl, sep='\t', usecols=['reads_order', anno_col])
    df[fid] = df.apply(lambda x: 2 if 'JRM' in x['reads_order'] else 1, axis=1)
    return df.groupby(anno_col).sum()[[fid]].reset_index().sort_values(fid, ascending=False)


def split_KOs(tbl):
    old_list = tbl.values.tolist()
    
    new_list = []
    for item in old_list:
        if ';' in item[0]:
            for KO in item[0].split(';'):
                new_list.append([KO, item[1]])
        else:
            new_list.append(item)
    df = pd.DataFrame(new_list, columns=tbl.columns)
    return df.groupby(tbl.columns[0]).sum().reset_index().sort_values(tbl.columns[1], ascending=False)


if __name__ == '__main__':
    arg_parser = ArgumentParser(
        description='Extract kegg abundance from IGC annotation of samples')

    arg_parser.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='IGC annotation of samples')
    arg_parser.add_argument(
        '-f',
        '--fid',
        dest='fid',
        required=True,
        help='file id')
    arg_parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output file')

    args = arg_parser.parse_args()

    tbl = args.input
    fid = args.fid
    output = args.output
    
    df = kegg_output(tbl, fid)
    df_new = split_KOs(df)
    df_new.to_csv(output, index=False, sep='\t')
