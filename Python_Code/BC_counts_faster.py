#!/usr/bin/python
# coding: utf-8


##### Barcode Counts (Faster)
# This still needs work but avoids looping through the dictionary
# Finds Barcode based on position

import pandas as pd

Barcodes_df = pd.read_csv("/active/cherry_t/091923_LVBM_data/resources/Twist_barcodes_RC.txt", names=['Barcode', 'ID'], sep="\t")

def BC_counts_position(path, barcodes):# internally defined dictionary
    Barcodes_dict_count = dict.fromkeys(barcodes, 0)
    line_no = 0
    with open(path, "rt") as fq:# loop through the fastq
        for line in fq:
            if line_no % 4 == 1:# only pull sequence
                seq = line[28:46]# anticipated BC position
                if seq in Barcodes_dict_count:# just look to see if the sequence at position 'seq' is a known barcode
                    Barcodes_dict_count[seq] +=1
            line_no +=1
    return(Barcodes_dict_count)

# eventually I can take these paths from the command line and loop through them
counts_1 = BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/1.Twist_cDNA_S1_L001_R1_001.fastq", Barcodes_df['Barcode'])
counts_2 = BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/1.Twist_cDNA_S1_L001_R2_001.fastq", Barcodes_df['Barcode'])
counts_3 = BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/2.Twist_cDNA_S1_L002_R1_001.fastq", Barcodes_df['Barcode'])
counts_4 = BC_counts_position("/active/cherry_t/091923_LVBM_data/Twist_cDNA/2.Twist_cDNA_S1_L002_R2_001.fastq", Barcodes_df['Barcode'])

Barcodes_df['Counts_1'] = Barcodes_df['Barcode'].map(counts_1)
Barcodes_df['Counts_2'] = Barcodes_df['Barcode'].map(counts_2)
Barcodes_df['Counts_3'] = Barcodes_df['Barcode'].map(counts_3)
Barcodes_df['Counts_4'] = Barcodes_df['Barcode'].map(counts_4)
Barcodes_df['Combined_Counts'] = Barcodes_df['Counts_1'] + Barcodes_df['Counts_2'] + Barcodes_df['Counts_3'] + Barcodes_df['Counts_4']

Barcodes_df.to_csv("/active/cherry_t/091923_LVBM_data/Twist_cDNA/Twist_cDNA_Barcode_counts_alt.txt", sep='\t', index=False)

