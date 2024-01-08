# This script is producing the "dirty" ML tree.

library(BiocManager)
library(BiocGenerics)
library(msa)
library(seqinr)
library(ape) 
library(Biostrings)
library(ade4)
library(readr)
library(microseq)
library(phangorn)

# First of all, we should get the sequences to prepare the *.fasta file:

start_csv <- read_csv("start_csv.csv")
numbers <- as.data.frame(t(start_csv), col.names = c('0'))
sequences <- read.GenBank(numbers)
attr(sequences, "species")                     
str(sequences)
names <- paste(names(sequences), attr(sequences, "species"), sep="_")
write.dna(sequences, file="sequences_fst.fasta", format = "fasta")
sequences_seqinr <- read.fasta(file="sequences_fst.fasta", seqtype = "DNA", as.string = TRUE)
write.fasta(sequences = sequences_seqinr, names = names, file.out = "final_output.fasta")

fn <- "sequences_fst.fasta"
 if (file.exists(fn)) {
   file.remove(fn)
 }

 # During the next step, we should prepare the MSA:

Sequences_for_alignment <- readDNAStringSet("final_output.fasta")
myMuscleAlignment <- msa(Sequences_for_alignment, "Muscle", type = "dna")
writeXStringSet(as(unmasked(myMuscleAlignment), "XStringSet"), file="./myMuscleAlignment.fasta")

# And now we can trim them:
Sequences_for_trimming <- readFasta("myMuscleAlignment.fasta")
trimmedMuscleAlignment <- msaTrim(Sequences_for_trimming, gap.end = 0.5, gap.mid = 0.9)
writeFasta(trimmedMuscleAlignment, "trimmedMuscleAlignment.fasta", width = 0)

# Last step - tree building:
fasta <- read.phyDat("trimmedMuscleAlignment.fasta", format = 'fasta')
#mt <- modelTest(fasta)
#best_model <- as.pml(mt)
#best_model
tree <- pml_bb(fasta, model="GTR+G(4)+I")
write.tree(tree$tree, "ML_tree.tre")