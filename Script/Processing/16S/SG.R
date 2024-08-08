## Script run on Apphub account of uni Greifswald
## Processing of 'Salinity_Gradient'

# Loading necessary libraries
library(dada2)
library(DECIPHER)
library(phyloseq)
library(speedyseq)
library(ggplot2)
library(Biostrings)
library(ape)
library(microbiome)
library(vegan)

set.seed(68)

# Load the merged sequence table for further processing
seqtab_Sal <- readRDS("/home/rstudio/Documents/Seqtab/Seqtab_16S_Salinity.rds")

# Display the dimensions of the sequence table
dim(seqtab_Sal)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab_Sal)))

# Our expected amplicon size is 250 bp, so we discard all samples that are too short or too long, allowing for a 10 bp wiggle room either way
# The amplified V4 region of the 16S rRNA gene does not have great length variability
seqtab_cut_Sal <- seqtab_Sal[,nchar(colnames(seqtab_Sal)) %in% 240:260]

# Verify the distribution of sequence lengths after trimming
table(nchar(getSequences(seqtab_cut_Sal)))

# Collapse no mismatch sequences with a minimum overlap of 200 bp
seqtab_cut_uni_Sal <- collapseNoMismatch(seqtab_cut_Sal, minOverlap = 200, verbose = TRUE)

# Remove chimeric sequences
seqtab.nochim_Sal <- removeBimeraDenovo(seqtab_cut_uni_Sal, method="consensus", multithread=TRUE, verbose=TRUE)

# Create a DNAStringSet from the ASVs
dna_Sal <- DNAStringSet(getSequences(seqtab.nochim_Sal))  
names(dna_Sal) <- DNAStringSet(dna_Sal)

# Create a phyloseq object from the sequence (ASV)-abundance-table and the DNA stringset
phy_16S_Sal <- phyloseq(otu_table(seqtab.nochim_Sal, taxa_are_rows=FALSE), refseq(dna_Sal))

# Change ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(phy_16S_Sal) <- paste0("ASV", seq(ntaxa(phy_16S_Sal)))

# Export ASV and sequence table for classification with CREST4
writeXStringSet(refseq(phy_16S_Sal), "Seqs_16S_Sal.fasta", format = "fasta")
asv_table_df_t_Sal <- t(as.data.frame(otu_table(phy_16S_Sal)))
write.csv(asv_table_df_t_Sal, "Asv_16S_counts_Sal.csv")

# Save the workspace
save.image("aphy_16S_Sal.RData")

# After running CREST4 on Jupyter for classification

# Define taxonomic ranks
tax_ranks <- c("ASV", "root", "type", "kingdom", "superkingdom", "superphylum", "phylum", "class", "order", "family", "genus", "species")

# Read and clean the classification results
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/Documents/fasta_crest4/all_seqs_16S_Sal.fasta.crest4/assignments_16S.txt")), 
                     header = FALSE, fill = TRUE, na.strings = "", sep = ";", col.names = tax_ranks)

# Fill in NA columns with the entry from the previous cell (last classified taxonomic rank) and add "unclassified"
for (i in 1:ncol(assign)) {
  assign[,i] <- ifelse(is.na(assign[,i]), paste0(assign[,i-1], "_unclassified"), assign[,i])
}

# Set row names to ASV names
row.names(assign) <- assign$ASV

# Convert the assignment data frame to a matrix
assign_ma <- as.matrix(assign[,-1])
assign_ma <- gsub("(_unclassified_.{1,})", "_unclassified", assign_ma) # Remove multiple occurrences of "_unclassified"

# Create a taxonomic table
tax_phy <- tax_table(assign_ma)

# Read in the OTU table
pro_otu <- read.csv("~/Documents/counts/Asv_16S_counts_Sal.csv", header = TRUE, row.names = 1, sep= ",")

# Transpose the OTU table to put sample names in rows
pro_otu <- t(pro_otu)

# Extract sample names
sample_names <- rownames(pro_otu)

# Remove the prefix "X" from sample names
sample_names <- sub("^X", "", sample_names)

# Replace dots with hyphens in sample names
sample_names <- gsub("\\.", "-", sample_names)

# Extract short sample names
sample_names <- sapply(strsplit(sample_names, "515F-806R-EMP2016-P03-"), `[`, 2)
short_names <- as.vector(t(sample_names))
row.names(pro_otu) <- short_names
rm(sample_names)

# Read in metadata
meta <- read.csv("~/nextcloud/New_Meta_Sal.csv", header = TRUE, sep = ",", na.strings = c(""," ","NA"), nrows = 26, as.is = TRUE) 
row.names(meta) <- meta$SN

# Create an OTU table
pro_OTU <- otu_table(pro_otu, taxa_are_rows=FALSE)

# Merge phyloseq objects
phys_Sal <- merge_phyloseq(pro_OTU, tax_phy, phy_16S_Sal, sample_data(meta))
sample_names(phys_Sal) <- sample_data(phys_Sal)$Sample_ID
print(phys_Sal)

# Save the workspace
save.image("~/Documents/workspace/16S_AM_Sal")

# Define the output directory
output_dir <- "/home/rstudio/Documents/phyloseq"

# Save the phyloseq object in RDS format
saveRDS(phys_Sal, file = paste0(output_dir, "/Phyloseq_Sal_16S.rds"))