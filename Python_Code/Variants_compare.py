#!/usr/bin/env python3

"""
Variants Comparison Pipeline
This is a series of functions to compare expression of barcodes from the wt sequences to differing levels of mutational load
Additionally, it outputs fasta files for running dsvm
Eventually I will get this stript to run dsvm and append the results to the individual comparisons table

Usage: Just run as is. I could make this with command line inputs for the files, but it will be a lot of variables.
For now it's hard coded, until I make it prettier.
"""

import os
import re
import pandas as pd
import pybedtools
from scipy import stats
import subprocess

# Output basic summary of number of barcodes per enhancer('ID')
def barcodes_summary(tables, table_names, output):
    summary_tables = {}
    # You could do this with a list of paths
    # But I wanted more control over the names of the categories
    for i in range(len(tables)):
        table = tables[i]
        counts = table['ID'].value_counts()
        df_name = table_names[i]
        summary_tables[df_name] = counts
    summary_df = pd.concat(summary_tables, axis=1)
    summary_df = summary_df.fillna(0)
    summary_df.to_csv(output, index=True, sep='\t')

# Output the average change in expression from normalized counts between mutational categories
def variants_comparisons_bulk(novariants_filtered, mutagenized_dfs, mutagenized_names, output):
    comparison_results = pd.DataFrame(columns=['ID', 'Mutagenized_Category', 'Change', 'P_Value'])
    for id in novariants_filtered['ID'].unique():
# filter WT expression table by enhancer
        novariants_id = novariants_filtered[novariants_filtered['ID'] == id]
        # Same overcomplicated laziness
        for i in range(len(mutagenized_dfs)):
            df = mutagenized_dfs[i]
            category = mutagenized_names[i]
            mut_id = df[(df['ID'] == id)] # filter again
            # get stats and change of expression
            t_statistic, p_value = stats.ttest_ind(df['Normalized_Counts'], novariants_id['Normalized_Counts'])
            mean_change = df['Normalized_Counts'].mean() - novariants_id['Normalized_Counts'].mean()
            new_row = pd.DataFrame({'ID':id, 'Mutagenized_Category':category, 'Change': mean_change, 'P_Value': p_value}, index=[0])
            comparison_results = pd.concat([comparison_results, new_row], ignore_index=True)
    comparison_results.to_csv(output, sep='\t', index = False)

# Because each variant will have a different impact
# Get changes per barcode in mutants
def variants_comparisons_all(novariants_filtered, mutagenized_dfs, mutagenized_names):
    individual_comparisons = pd.DataFrame(columns=['Barcode', 'Variant_Calls', 'ID', 'Mutagenized_Category', 'Normalized_Counts', 'Xn_Change'])
    for id in novariants_filtered['ID'].unique():
        novariants_id = novariants_filtered[novariants_filtered['ID'] == id]
        mean_xn = novariants_id['Normalized_Counts'].mean()
        for i in range(len(mutagenized_dfs)):
            df = mutagenized_dfs[i]
            cat = mutagenized_names[i]
            mut_id = df[(df['ID'] == id)]
            for index, line in mut_id.iterrows():
                Xn_Change = line['Normalized_Counts'] - mean_xn
                new_row = pd.DataFrame({'Barcode':line['Barcode'], 'Variant_Calls': line['Variant_Calls'], 'ID': id, 'Mutagenized_Category': cat, 'Normalized_Counts': line['Normalized_Counts'], 'Xn_Change': Xn_Change}, index=[0])
                individual_comparisons = pd.concat([individual_comparisons, new_row], ignore_index=True)
    #individual_comparisons.to_csv(output, sep='\t', index=False)
    return individual_comparisons

# This function could be wrapped into the next but I was getting annoyed with troubleshooting
def get_fasta(bed, ref_fa='/active/cherry_t/Leah/Resources/SEQ_UTILS/hg38.fa'):
    position_seq = bed.sequence(fi=ref_fa)
    fasta_str = open(position_seq.seqfn).read()
    # This should be modified if you are expecting more than one sequence
    fasta_list = fasta_str.split('\n')[1:-1] # skip header and empty row
    return ''.join(fasta_list)

# From a list of sequences and snps in each, output fastas designed for dSVM
# This is designed to take output from variants_comparisons_all
def variant_fasta(table, ref_out, var_out, ref_fa='/active/cherry_t/Leah/Resources/SEQ_UTILS/hg38.fa'):
    # Because opening in append mode, clear any preexisting files
    if os.path.exists(ref_out):
        os.remove(ref_out)
    if os.path.exists(var_out):
        os.remove(var_out)
    num_mismatches = 0
    num_correct = 0
    # iterate through individual comparisons table
    with open(ref_out, 'a') as wt_out, open(var_out, 'a') as mut_out:
        
        for index, line in table.iterrows():
            barcode = line['Barcode'] # reserve BC for sequence names
            # Grab original position from enh name, get seq
            chrom, start, stop = re.search(r'(chr\w+):(\d+)-(\d+)', line['ID']).groups()
            if chrom and start and stop:
                bed = pybedtools.BedTool(f'{chrom}\t{start}\t{stop}', from_string=True)
                sequence_string = get_fasta(bed, ref_fa).upper()
                variants = line['Variant_Calls'].split(';')
                variant_data = [re.search(r'(\d+)(\w)>(\w)', x).groups() for x in variants]
                # make substitutions along the sequence according to variants per barcode
                alt_str = sequence_string
                for pos, ref, alt in variant_data:
                    str_pos = int(pos) - int(start)
                    if str_pos < 0 or str_pos > 245:
                        #print(f'Position {str_pos} out of range, skipping')
                        num_mismatches +=1
                    elif sequence_string[str_pos-1] == ref:
                        alt_str = alt_str[:str_pos-1] + alt + alt_str[str_pos:]
                        num_correct += 1
                    else:
                        #print(f'Mismatch for position {pos} in {chrom}:{start}-{stop}. string pos = {str_pos}. Ref: {ref}, Alt: {alt}, found nt: {sequence_string[str_pos-1]}')
                        num_mismatches += 1
                        #raise ValueError("Reference nucleotide doesn't match, fix your function")
                # build new entry in fasta format
                out_ref = f'>{barcode}\n{sequence_string}\n'
                out_alt = f'>{barcode}\n{alt_str}\n'
                wt_out.write(out_ref)
                mut_out.write(out_alt)
            else:
                raise ValueError('ID format does not contain positional information')
    print(f'Correct substitutions made: {num_correct}\nIncorrect substitutions skipped: {num_mismatches}')

# Next function
# take prior fastas, model of interest
# run dsvm (to output table)
def run_dsvm(model, ref_fa, var_fa, output):
    command = f'perl /active/cherry_t/Leah/Resources/deltasvm_script/deltasvm.pl {ref_fa} {var_fa} {model} {output}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(result.stdout)
    else:
        print("Error:", result.stderr)

def get_dsvm_dict(dsvm_out):
    dsvm_dict = {}
    with open(dsvm_out, 'r') as dsvm:
        for line in dsvm:
            key, value = line.strip().split('\t')
            dsvm_dict[key] = value
    return dsvm_dict

def main():
    
    # Ideally this would take things from the command line
    # but there's a lot of distinct files here so
    # hard code it, sorry

    # All reads from BC Analysis split by mutational load
    novariants_raw = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_novariant_pysam_BC_table_raw.txt", sep='\t')
    onevariant_raw = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_1variant_pysam_BC_table_raw.txt", sep='\t')
    fewvariants_raw = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_fewvariants_pysam_BC_table_raw.txt", sep='\t')
    manyvariants_raw = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_manyvariants_pysam_BC_table_raw.txt", sep='\t')
    empty_raw = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_empty_pysam_BC_table_raw.txt", sep='\t')
    
    # Reads filtered by pDNA >= 10 and normalized
    novariants_filtered = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_novariant_pysam_BC_table_above10.txt", sep='\t')
    onevariant_filtered = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_1variant_pysam_BC_table_above10.txt", sep='\t')
    fewvariants_filtered = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_fewvariants_pysam_BC_table_above10.txt", sep='\t')
    manyvariants_filtered = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_noindel_manyvariants_pysam_BC_table_above10.txt", sep='\t')
    empty_filtered = pd.read_csv("/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_empty_pysam_BC_table_above10.txt", sep='\t')
    
    # Tables and names for iteration
    tables = [novariants_raw, novariants_filtered, onevariant_raw, onevariant_filtered, fewvariants_raw, fewvariants_filtered, manyvariants_raw, manyvariants_filtered, empty_raw, empty_filtered]
    table_names =  ['novariants_raw', 'novariants_filtered', 'onevariant_raw', 'onevariant_filtered', 'fewvariants_raw', 'fewvariants_filtered', 'manyvariants_raw', 'manyvariants_filtered', 'empty_raw', 'empty_filtered']
    
    mutagenized_dfs = [onevariant_filtered, fewvariants_filtered, manyvariants_filtered]
    mutagenized_names = ['onevariant_filtered', 'fewvariants_filtered', 'manyvariants_filtered']
    
    # Model for running
    # Which is actually 11mer scores from a model
    bulk_model = '/active/cherry_t/Leah/Analyses/gkm_svm/RandNeg_analysis/Kernel2/11mers/all_11mers_merged_ATAC_ChIP_randneg_param_gkmpredict.txt'
    rods_model = '/active/cherry_t/Leah/Analyses/gkm_svm/scATAC/11mers/all_11mers_Mature_Rods_noMatIn_t25k_.txt'
    mg_model = '/active/cherry_t/Leah/Analyses/gkm_svm/scATAC/11mers/all_11mers_Mature_Mullers_noMatIn_t25k_.txt'

    # Output names
    output_motif = '/active/cherry_t/Leah/MPRA/Twist_WT_enhs/variants/pipeline_outs/Twist_variants_comparison'
    
    BC_counts_out = output_motif + 'BC_summary_table.txt'
    mut_comparison_bulk_out = output_motif + 'mutant_bulk_comparison_table.txt'
    mut_comparison_individual_out = output_motif + 'all_variants_comparison_table.txt'
    ref_out = output_motif + 'variants_wt.fa'
    mut_out = output_motif + 'variants_mut.fa'
    bulk_dsvm_out = output_motif + 'bulk_dsvm.txt'
    rod_dsvm_out = output_motif + 'rod_dsvm.txt'
    mg_dsvm_out = output_motif + 'mg_dsvm.txt'


    # Run functions
    barcodes_summary(tables, table_names, BC_counts_out)
    variants_comparisons_bulk(novariants_filtered, mutagenized_dfs, mutagenized_names, mut_comparison_bulk_out)
    individual_comparisons = variants_comparisons_all(novariants_filtered, mutagenized_dfs, mutagenized_names)
    variant_fasta(individual_comparisons, ref_out, mut_out)
    # run dsvm through models
    run_dsvm(bulk_model, ref_out, mut_out, bulk_dsvm_out)
    run_dsvm(rods_model, ref_out, mut_out, rod_dsvm_out)
    run_dsvm(mg_model, ref_out, mut_out, mg_dsvm_out)

    # map dsvm to individual comparisons table. Maybe change variants_comparisons_all and variant_fasta to deal with pandas df rather than this nonsense
    bulk_dsvm = get_dsvm_dict(bulk_dsvm_out)
    rod_dsvm = get_dsvm_dict(rod_dsvm_out)
    mg_dsvm = get_dsvm_dict(mg_dsvm_out)
    
    individual_comparisons['Bulk_dSVM'] = individual_comparisons['Barcode'].map(bulk_dsvm)
    individual_comparisons['Rod_dSVM'] = individual_comparisons['Barcode'].map(rod_dsvm)
    individual_comparisons['MG_dSVM'] = individual_comparisons['Barcode'].map(mg_dsvm)
    individual_comparisons.to_csv(mut_comparison_individual_out, sep='\t', index=False)

if __name__ == "__main__":
    main()


