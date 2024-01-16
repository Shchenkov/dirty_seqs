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
library(catGenes)

# First of all, we should get the sequences to prepare the *.fasta file:
start_csv <- read_csv("start_csv.csv")
numbers <- as.data.frame(t(start_csv), col.names = c('0'))
sequences <- read.GenBank(numbers)
attr(sequences, "species")                     
str(sequences)
names <- paste(names(sequences), attr(sequences, "species"), sep="_")
names <- gsub("-", "_", names)
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

# And now we can trim the msa:
Sequences_for_trimming <- readFasta("myMuscleAlignment.fasta")
trimmedMuscleAlignment <- msaTrim(Sequences_for_trimming, gap.end = 0.5, gap.mid = 0.9)
writeFasta(trimmedMuscleAlignment, "trimmedMuscleAlignment.fasta", width = 0)

# To save trimmed alignment as *.nexus:
fasta_for_nexus <- readFasta('trimmedMuscleAlignment.fasta')
fasta_for_nexus_df <- as.data.frame(fasta_for_nexus)
nexusdframe(fasta_for_nexus_df, "trimmedMuscleAlignment.nxs")

# To test the best fitted model:
mt <- modelTest(fasta)
best_model <- as.pml(mt)
best_model
mt_df <- as.data.frame(mt)
write.csv(mt_df, file = "model.testing.csv")

# And to write Bayes-block to the end of the nexus-file:
bayes_block = "begin mrbayes;
set autoclose=yes nowarn=yes;
lset covarion=no nucmodel=4by4 nst=6 ngammacat=8 rates=Invgamma; mcmcp ngen=15000000 nruns=2 nchains=4 samplefreq=100 printf=1000 stoprule=n stopval=0.01;
mcmc;
END;"
write(bayes_block, file = "trimmedMuscleAlignment.nxs", append = TRUE)

# Last step - tree building:
fasta <- read.phyDat("trimmedMuscleAlignment.fasta", format = 'fasta')
tree <- pml_bb(fasta, model="GTR+G(4)+I")
write.tree(tree$tree, "ML_tree.tre")
