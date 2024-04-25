#!/usr/bin/env python3

"""
Barcode Analysis Pipeline
Counts defined barcodes in cDNA and pDNA fastq files to output counts tables and barplot.

Usage:
    BC_analysis.py --cdna <cDNA_INPUT> --pdna <pDNA_INPUT> --barcodes <BARCODES> --out <OUTPUT_PREFIX>

Options:
    --cdna <cDNA_INPUT>         Input cDNA fastq.gz file(s)
    --pdna <pDNA_INPUT>         Input pDNA fastq.gz file(s)
    --barcodes <BARCODES>       Barcode table (Be aware of reverse complement, format)
    --out <OUTPUT_PREFIX>       Prefix for all output files
"""

import argparse
from Bio import SeqIO
import sys
import re
import pandas as pd
import gzip
import seaborn as sns
import matplotlib.pyplot as plt


def BC_counts(path, counts_dict):
    with gzip.open(path, 'rt') as fq:
        for record in SeqIO.parse(fq, 'fastq'):
            read = str(record.seq).upper()
            # Find BC by adjacent seq
            # Edit adjacent seq as needed
            # Edit BC RegEx as needed
            findBC = re.search(r'GGCAGAGGGAAAAAGATC(([AT][GC]){9})', read)
            if findBC:
                BC = findBC.group(1)
                if BC in counts_dict:
                    counts_dict[BC] += 1
    return counts_dict

def plot_counts(table_above10, out_plot, grouping = 'median'):
    table = table_above10.copy()
    table[grouping] = table.groupby('ID')['Normalized_Counts'].transform(grouping)
    table.sort_values(grouping, ascending=False, inplace=True)
    sns.set(style='ticks')
    plt.figure(figsize=(12,6))
    plot_all = sns.boxplot(x='ID', y='Normalized_Counts', data = table)
    plt.yscale('log')
    plot_all.tick_params(axis = 'x', rotation=90)
    plt.xlabel('Enhancer')
    plt.ylabel('Normalized Counts (cDNA/pDNA)')
    plt.savefig(out_plot, format = 'pdf', bbox_inches = 'tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Counts Barcodes in cDNA and pDNA and outputs a normalized table and plot')
    parser.add_argument('--cdna', metavar='cDNA_INPUT', type=str, nargs='+', help='Input cDNA fastq.gz files')
    parser.add_argument('--pdna', metavar='pDNA_INPUT', type=str, nargs='+', help='Input pDNA fastq.gz files')
    parser.add_argument('--barcodes', metavar='BARCODES', type=str, help='Barcode table (please use Reverse Complement)')
    parser.add_argument('--out', metavar='OUTPUT_PREFIX', type=str, help='Output prefix for output files')
    args = parser.parse_args()

    if args.barcodes:
        # Make master dictionary and table
        Barcodes_dict_count = {}
        with open(args.barcodes, 'r') as bc:
            table = pd.read_csv(bc, comment='#', sep='\t', header = None)
            table.columns = ['Barcode', 'Variant_Calls', 'ID']
    Barcodes_dict = {barcode: 0 for barcode in table['Barcode']}

    cDNA_Barcodes = Barcodes_dict.copy()
    if args.cdna:
        for fq in args.cdna:
            # Process barcodes to main counts dictionary
            cDNA_Barcodes = BC_counts(fq, cDNA_Barcodes)

    pDNA_Barcodes = Barcodes_dict.copy()
    if args.pdna:
        for fq in args.pdna:
            # Process barcodes to main counts dictionary
            pDNA_Barcodes = BC_counts(fq, pDNA_Barcodes)

    if args.out:
        # make output file names
        out_table = str(f'{args.out}_BC_table_raw.txt')
        out_table_above10 = str(f'{args.out}_BC_table_above10.txt')
        out_plot = str(f'{args.out}_above10_boxplot.pdf')

    # map count dictionaries
    table['cDNA_counts'] =  table['Barcode'].map(cDNA_Barcodes)
    table['pDNA_counts'] = table['Barcode'].map(pDNA_Barcodes)

    # output raw counts
    table.to_csv(out_table, sep="\t", index=False)

    # filter by pDNA counts
    Counts_above10 = table.loc[table['pDNA_counts'] >= 10].copy()
    # normalize counths
    Counts_above10['Normalized_Counts'] = Counts_above10['cDNA_counts'] / Counts_above10['pDNA_counts']
    # output normalized counts
    Counts_above10.to_csv(out_table_above10, sep="\t", index=False)
    # plot
    plot_counts(Counts_above10, out_plot)

if __name__ == "__main__":
    main()
