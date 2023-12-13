# importing windows:
with open("/home/lvand1/training_data/20p_outgroup_CRX_ALL_crossref.bed", "r") as fi:
    windows = fi.read().splitlines()

# breaking up windows into 1bp coords in bed format:
nt_bed = []
for i in range(len(windows)):
    bed = windows[i].split("\t")
    for j in range(int(bed[1]),int(bed[2])):
        nt_bed.append(bed[0] + '\t' + str(j) + '\t' + str(j+1))

# write to output file:
fo = open("/home/lvand1/training_data/CRX20p_1bpsplit.bed", "w")
fo.write("\n".join(nt_bed))
fo.close()
