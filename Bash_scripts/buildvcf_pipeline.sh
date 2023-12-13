#!/bin/bash

#Processing merged dSVM bedgraph to create VCF
# Input 1: bedGraph
# Input 2: Output file

# Make tmp file:
cat $1 > tmp_cat.bedGraph

#Save IDs
awk -v OFS="\t" '{print $4}' tmp_cat.bedGraph > tmp_dSVM_IDs.txt

#Isolate Ref and Alt BPs
sed -i 's/:/\t/g' tmp_cat.bedGraph 
sed -i 's/->/\t/g' tmp_cat.bedGraph 
sed -i 's/--/\t/g' tmp_cat.bedGraph

#Make column temp files
awk -v OFS="\t" '{print $1, $2}' tmp_cat.bedGraph > tmp_dsvm_trim.txt
awk -v OFS="\t" '{print $6, $7}' tmp_cat.bedGraph > tmp_dsvm_vcf_cols.txt
#add extra vcf cols. Could be done in prior step I think.
awk -F '\t' '{$(NF+1)="." FS ".";}1' OFS='\t' tmp_dsvm_vcf_cols.txt > tmp_dsvm_vcf_cols_xtra.txt
#store scores
awk -v OFS="\t" '{print $9}' tmp_cat.bedGraph > tmp_dsvm_scores_only.txt
sed -i -e 's/^/DSVM=/' tmp_dsvm_scores_only.txt

#Combine
#If IDs not desired, in extra step, use awk to add dots into third column of "trim" file
paste tmp_dsvm_trim.txt tmp_dSVM_IDs.txt tmp_dsvm_vcf_cols_xtra.txt tmp_dsvm_scores_only.txt > tmp_dsvm_vcfish.txt

#sort it.
module load bedtools/2.27.1
bedtools sort -i tmp_dsvm_vcfish.txt > tmp_sorted_dsvm_vcfish.txt

#Add a header
cat ~/deltasvm/vcf/vcf_header.txt tmp_sorted_dsvm_vcfish.txt > tmp_sorted_dsvm.vcf

#Zip & index
module load BCFtools/1.9-foss-2018b
bgzip -c tmp_sorted_dsvm.vcf > $2
tabix -p vcf $2

#rm *tmp*

exit
