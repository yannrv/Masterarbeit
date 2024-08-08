#### Script run on Apphub account of uni Greifswald
## Pre-processing of 'Seastore_Project'

# Loading necessary libraries
library(dada2)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(phangorn)

set.seed(68)

# Define the path to the sequence data
path <- "/home/rstudio/Downloads/my data seq/dada2_in_18S"

# List all files in the specified path
list.files(path)

# Sort files to ensure forward/reverse reads are in the same order
fnFs <- sort(list.files(path, pattern ="_R1.fastq"))
fnRs <- sort(list.files(path, pattern ="_R2.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Specify the full path to the forward and reverse reads
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Visualizing the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnFs[35:51])
plotQualityProfile(fnFs[54:70])

# Quality scores for specific forward reads
qualscoresF <- plotQualityProfile(fnFs[c(6,30,54,55)])
qualscoresF + scale_x_continuous(breaks=seq(0,300,25))

# Visualizing the quality profiles of the reverse reads
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs[34:51])
plotQualityProfile(fnRs[54:70])

# Quality scores for specific reverse reads
qualscoresR <- plotQualityProfile(fnRs[c(6,30,52,55)])
qualscoresR + scale_x_continuous(breaks=seq(0,300,25))

# Create another path for filtered files
path2 <- "/home/rstudio/Downloads"

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path2, "filtrer_18S_Seastore", paste0(sample.names, "_F_filt.fastq.bz2"))
filtRs <- file.path(path2, "filtrer_18S_Seastore", paste0(sample.names, "_R_filt.fastq.bz2"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim the sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=175, truncLen=c(210,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
head(out)

# Learn the error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# Visualize error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Perform sample inference for forward and reverse reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Dereplication of sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data frame from the first sample
head(mergers[[1]])

# Construct the sequence table
seqtab <- makeSequenceTable(mergers)

# Check the dimensions of the sequence table
dim(seqtab)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Define the output directory
output_dir <- "/home/rstudio/Documents/Seqtab"

# Save the sequence table in RDS format
saveRDS(seqtab, file = paste0(output_dir, "/Seqtab_18S_Seastore.rds"))