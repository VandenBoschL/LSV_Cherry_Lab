
##Write a reverse complement of a given sequence in a dataframe
#This assumes that the barcodes are in the first column of the table.
# Before running, in command line run "ml R/4.1.0-foss-2020b"
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_path <- args[1]

output <- args[2]

rev_comp <- function(seq){
  complement = c("A"="T", "T"="A", "G"="C", "C"="G")
  reverse_seq = as.character(toupper(rev(strsplit(seq, split="")[[1]])))
  revcomp_seq = complement[reverse_seq]
  return(paste(revcomp_seq, collapse=""))
}

dataframe <- read.delim(file_path, sep = "\t", header= FALSE)

dataframe$V3 <- sapply(dataframe$V1, rev_comp)

out_df <- data.frame(V1 = dataframe$V3, V2= dataframe$V2)
write.table(out_df, file= output, sep="\t", col.names = FALSE, row.names=FALSE, quote = FALSE)


