#!/usr/bin/python
# coding: utf-8

#### Barcode Counts
# This is the original script for finding MPRA barcodes in RNA-seq data
# Currently hard coded because I'm still developing this function
# Finds barcodes by position only
# Takes a loooong time to run if there's too many barcodes

import pandas as pd

Barcodes_df = pd.read_csv("/active/cherry_t/091923_LVBM_data/resources/Twist_barcodes_RC.txt", names=['Barcode', 'ID'], sep="\t")
Barcodes_dict_count = dict.fromkeys(Barcodes_df['Barcode'], 0)

def BC_counts_position(path):
    line_no = 0
    with open(path, "rt") as fq:
        for line in fq:# loop through the fastq
            if line_no % 4 == 1:# only pull sequence
                for barcode in Barcodes_dict_count.keys():# fix this to not have a hardcoded dictionary
                    if barcode == line[28:46]:# anticipated BC position
                        Barcodes_dict_count[barcode] +=1
            line_no +=1

# Because this is so hard coded, every iteration of the function appends counts to the dictionary
BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/1.Twist_cDNA_S1_L001_R1_001.fastq")
BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/1.Twist_cDNA_S1_L001_R2_001.fastq")
BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/2.Twist_cDNA_S1_L002_R1_001.fastq")
BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/2.Twist_cDNA_S1_L002_R2_001.fastq")

Barcodes_df['Count'] = Barcodes_df['Barcode'].map(Barcodes_dict_count)

Barcodes_df.to_csv("/active/cherry_t/091923_LVBM_data/Twist_cDNA/Twist_cDNA_Barcode_counts.txt", sep='\t', index=False)

