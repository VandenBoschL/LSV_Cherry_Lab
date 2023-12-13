
##Write a reverse complement of a given sequence in a dataframe
#This does not assume a position for the sequences and will return a df with the sequences replaced
# Before running, in command line run "ml R/4.1.0-foss-2020b"
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_path <- args[1]
sequence_position <- args[2]
output <- args[3]

rev_comp <- function(seq){
  complement = c("A"="T", "T"="A", "G"="C", "C"="G")
  reverse_seq = as.character(toupper(rev(strsplit(seq, split="")[[1]])))
  revcomp_seq = complement[reverse_seq]
  return(paste(revcomp_seq, collapse=""))
}

dataframe_seqs <- read.delim(file_path, sep = "\t", header= FALSE)

loc <- paste0("V", sequence_position)

RC_sequences <- sapply(dataframe_seqs[, loc], rev_comp)

out_df <- dataframe_seqs
out_df[,loc] <- RC_sequences

write.table(out_df, file= output, sep="\t", col.names = FALSE, row.names=FALSE, quote = FALSE)

