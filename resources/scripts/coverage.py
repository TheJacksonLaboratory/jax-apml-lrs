"""
Collects mean sequencing depth per chromosome from samtools coverage output,
computes the autosomal average, and saves results as a TSV file.

Usage:
    python3 coverage.py --covfile /path/to/samtools_coverage.txt

Output:
    meandepth.tsv  -- two-column TSV with per-chromosome and average mean depth
"""

import argparse
import pandas as pd


def coverage(file):
    """
    Read samtools coverage output, extract autosomal mean depth,
    append genome-wide average row, and write to meandepth.tsv.
    """
    columns = [
        'chromosome', 'startpos', 'endpos', 'numreads',
        'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq'
    ]
    df = pd.read_csv(file, names=columns, sep=r'\s+')
    df = df.iloc[1:].reset_index(drop=True)
    df['meandepth'] = pd.to_numeric(df['meandepth'], errors='coerce')

    autosomes = ['chr' + str(x) for x in range(1, 23)]
    df = df[['chromosome', 'meandepth']]
    df_autosomes = df.loc[df['chromosome'].isin(autosomes)]

    avg_row = {'chromosome': 'average', 'meandepth': df_autosomes['meandepth'].mean()}
    out = pd.concat([df_autosomes, pd.DataFrame([avg_row])], ignore_index=True)
    out.to_csv('meandepth.tsv', sep='\t', index=False)


def get_args():
    parser = argparse.ArgumentParser(
        description='Compute mean sequencing depth from samtools coverage output.'
    )
    parser.add_argument(
        '--covfile',
        type=str,
        required=True,
        help='Path to samtools coverage output file'
    )
    return parser.parse_args()


def main():
    args = get_args()
    coverage(args.covfile)


if __name__ == '__main__':
    main()