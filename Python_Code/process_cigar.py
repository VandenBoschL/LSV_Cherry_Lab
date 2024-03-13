#!/usr/bin/env python3

##############################
### Process CIGAR and MD:Z ###
##############################

# Developed for Nextflow barcode pairing pipeline
# Per read in sam file, builds "human readable" info for variants
# Allows for downstream filtration of barcodes by sequence fidelity



import re

def process_cigar(cigar, start_pos):
    # split codes
    cigar_values = re.findall(r'(\\d+|[MDINSHP=X])',cigar)
    # include if statement to divert code list len == 2 to "no indels"
    if len(cigar_values) <= 2:
        return "No_Indels"
    # generate counter and output list
    else:
        readable_codes = []
        seq_length = 0
        for position, value in re.findall(r'(\\d+)([MIDNSHP=X])', cigar):
            # Matches, only advance counter
            if value == 'M':
                seq_length += int(position)
            # cigar codes, spit out position and code + bp changes
            else:
                chrom_pos = start_pos + seq_length
                readable_codes.append(f'{chrom_pos}{value}{position}bp')
        return ';'.join(readable_codes)

def process_mdz(mdz, start_pos, sequence):
    # split codes
    mdz_variants = re.findall(r'(\\d+)([A-Z]|\\^[A-Z]+)?', mdz)
    # include if statement to divert code list len == 1 to "no variants"
    if len(mdz_variants) == 1:
        return 'No_Variants'
    else:
        # generate counter and output list
        variants_output = []
        seq_length = 0
        deletions = 0
        # iterate through variants except last
        for position, variant in mdz_variants[:-1]:
            seq_position = int(position)
            seq_length += seq_position
            alt_base =  sequence[seq_length]
            genomic_pos = start_pos + seq_length + deletions
            seq_length += 1 # Account for variant
            if variant.startswith('^'):
                # Sequence deletion
                ref_value = variant[1:]
                variants_output.append(f'{genomic_pos}{ref_value}>*')
                deletions += len(ref_value) - 1 # Adjust future genomic_pos 
            else:
                variants_output.append(f'{genomic_pos}{variant}>{alt_base}')
        return ';'.join(variants_output)

