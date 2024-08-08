## Script run on Apphub account of Greifswald
## Processing of 'Hiddensea_Study'

# Load necessary libraries
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

# Load the sequence table
seqtab_HS <- readRDS("/home/rstudio/Documents/Seqtab/Seqtab_16S_HS.rds")

# Display the dimensions of the sequence table
dim(seqtab_HS)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab_HS)))

# Our expected amplicon size is 250 bp, so we discard all samples that are too short or too long, allowing for a 10 bp wiggle room either way
# The amplified V4 region of the 16S rRNA gene does not have great length variability
seqtab_cut_HS <- seqtab_HS[,nchar(colnames(seqtab_HS)) %in% 240:260]

# Verify the distribution of sequence lengths after trimming
table(nchar(getSequences(seqtab_cut_HS)))

# Collapse no mismatch sequences with a minimum overlap of 200 bp
seqtab_cut_uni_HS <- collapseNoMismatch(seqtab_cut_HS, minOverlap = 200, verbose = TRUE)

# Remove chimeric sequences
seqtab.nochim_HS <- removeBimeraDenovo(seqtab_cut_uni_HS, method="consensus", multithread=TRUE, verbose=TRUE)

# Create a DNAStringSet from the ASVs
dna_HS <- DNAStringSet(getSequences(seqtab.nochim_HS))  
names(dna_HS) <- DNAStringSet(dna_HS)

# Create a phyloseq object from the sequence (ASV)-abundance-table and the DNA stringset
phy_16S_HS <- phyloseq(otu_table(seqtab.nochim_HS, taxa_are_rows=FALSE), refseq(dna_HS))

# Change ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(phy_16S_HS) <- paste0("ASV", seq(ntaxa(phy_16S_HS)))

# Export ASV and sequence table for classification with CREST4
writeXStringSet(refseq(phy_16S_HS), "Seqs_16S_HS.fasta", format = "fasta")
asv_table_df_t_HS <- t(as.data.frame(otu_table(phy_16S_HS)))
write.csv(asv_table_df_t_HS, "Asv_16S_counts_HS.csv")

# Save the workspace
save.image("phy_16S_HS.RData")

## After running on Jupyter

# Define taxonomic ranks
tax_ranks <- c("ASV", "root", "type", "kingdom", "superkingdom", "superphylum", "phylum", "class", "order", "family", "genus", "species")

# Read and clean the classification results
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/Documents/fasta_crest4/Seqs_16S_HS.fasta.crest4/assignments.txt")), 
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

# Read in metadata
meta <- read.csv("~/nextcloud/New_Meta_HS.csv", header = TRUE, sep = ",", na.strings = c(""," ","NA"), nrows = 36, as.is = TRUE) 

# Set row names of metadata to sample IDs
row.names(meta) <- meta$Sample_ID

# Read in the OTU table
pro_otu <- read.csv("~/Documents/counts/all_asv_16S_counts_HS.csv", header = TRUE, row.names = 1, sep= ",")

# Transpose the OTU table
pro_otu <- t(pro_otu)

# Clean sample names by removing prefix "X" and replacing dots with hyphens
sample_names <- rownames(pro_otu)
sample_names <- sub("^X", "", sample_names)
sample_names <- gsub("\\.", "-", sample_names)
sample_names <- sapply(strsplit(sample_names, "515YF-B806R-"), `[`, 2)
row.names(pro_otu) <- sample_names
rm(sample_names)

# Create an OTU table
pro_OTU <- otu_table(pro_otu, taxa_are_rows=FALSE)

# Merge phyloseq objects
phys_HS <- merge_phyloseq(pro_OTU, tax_phy, phy_16S_HS , sample_data(meta))

# Print the phyloseq object
print(phys_HS)

# Save the workspace
save.image("~/Documents/workspace/16S_AM_HS")

# Define the output directory
output_dir <- "/home/rstudio/Documents/phyloseq"

# Save the phyloseq object in RDS format
saveRDS(phys_HS, file = paste0(output_dir, "/Phyloseq_HS_16S.rds"))
