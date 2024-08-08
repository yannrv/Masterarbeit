## Script run on Apphub account of uni Greifswald
## Processing of 'Seagrass_Vibrio'

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

# Merge run 18S sequence table
seqtab_V <- readRDS("/home/rstudio/Documents/Seqtab/Seqtab_18S_Vibrio.rds")

# Display the dimensions of the sequence table
dim(seqtab_V)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab_V)))

# Our expected amplicon size is 250 bp, so we discard all samples that are too short or too long, allowing for a 10 bp wiggle room either way
# The amplified V4 region of the 18S rRNA gene does not have great length variability
seqtab_cut_V <- seqtab_V[,nchar(colnames(seqtab_V)) %in% 240:260]

# Verify the distribution of sequence lengths after trimming
table(nchar(getSequences(seqtab_cut_V)))

# Collapse no mismatch sequences with a minimum overlap of 200 bp
seqtab_cut_uni_V <- collapseNoMismatch(seqtab_cut_V, minOverlap = 200, verbose = TRUE)

# Remove chimeric sequences
seqtab.nochim_V <- removeBimeraDenovo(seqtab_cut_uni_V, method="consensus", multithread=TRUE, verbose=TRUE)

# Create a DNAStringSet from the ASVs
dna_V <- DNAStringSet(getSequences(seqtab.nochim_V))  
names(dna_V) <- DNAStringSet(dna_V)

# Create a phyloseq object from the sequence (ASV)-abundance-table and the DNA stringset
phy_18S_Vi <- phyloseq(otu_table(seqtab.nochim_V, taxa_are_rows=FALSE), refseq(dna_V))

# Change ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(phy_18S_Vi) <- paste0("ASV", seq(ntaxa(phy_18S_Vi)))

# Export ASV and sequence table for classification with CREST4
writeXStringSet(refseq(phy_18S_Vi), "Seqs_18S_V.fasta", format = "fasta")
asv_table_df_t_V <- t(as.data.frame(otu_table(phy_18S_Vi)))
write.csv(asv_table_df_t_V, "Asv_18S_counts_V.csv")

# Save the workspace
save.image("phy_18S_Vi.RData")

## After running CREST4 on Jupyter for classification

# Read in the OTU table
pro_otu <- read.csv("~/Asv_18S_counts_V.csv", header = TRUE, row.names = 1, sep= ",")

# Transpose the OTU table to put sample names in rows
pro_otu <- t(pro_otu)

# Extract sample names
sample_names <- rownames(pro_otu)

# Remove the prefix "X" from sample names
sample_names <- sub("^X", "", sample_names)

# Replace dots with hyphens in sample names
sample_names <- gsub("\\.", "-", sample_names)

# Extract short sample names
sample_names <- sapply(strsplit(sample_names, "1183F-1443R-"), `[`, 2)
short_names <- as.vector(t(sample_names))
row.names(pro_otu) <- short_names
rm(sample_names)

# Correct specific sample names
rownames(pro_otu) <- gsub("GLO-SGY-2", "GLO-SGO-4_1:2", rownames(pro_otu))
rownames(pro_otu) <- gsub("GLO-SS", "GH7-SS", rownames(pro_otu))
rownames(pro_otu) <- gsub("SHV-SGY-2", "GLO-SGO-4_1:10", rownames(pro_otu))
rownames(pro_otu) <- gsub("SHV-WF", "GLO-SGO-4_1:20", rownames(pro_otu))
rownames(pro_otu) <- gsub("GH5-SGO", "GLO-SGO-4_1:100", rownames(pro_otu))

# Print the processed OTU table
print(pro_otu)

# Create an OTU table
pro_OTU <- otu_table(pro_otu, taxa_are_rows=FALSE)

# Define taxonomic ranks
tax_ranks <- c("ASV", "root", "type", "domain", "rank_I", "rank_II", "rank_III", "rank_IV", "rank_V", "rank_VI", "rank_VII", "rank_VIII", "rank_IX", "rank_X", "rank_XI")

# Read and clean the classification results
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/all_seqs_18S_V.fasta.crest4/assignments.txt")), 
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
meta <- read.csv("~/nextcloud/New_Meta_Vibrio.csv", header = TRUE, sep = ",", na.strings = c(""," ","NA"), nrows = 51, as.is = TRUE) 

# Set row names of metadata to sample IDs
row.names(meta) <- meta$Sample_ID

# Merge phyloseq objects
all_18S_phy_Vibrio <- merge_phyloseq(pro_OTU, tax_phy, sample_data(meta), phy_18S_Vi)

# Print the merged phyloseq object
print(all_18S_phy_Vibrio)

# Save the workspace
save.image("~/Documents/workspace/18S_AM_Vi")

# Define the output directory
output_dir <- "/home/rstudio/Documents/phyloseq"

# Save the phyloseq object in RDS format
saveRDS(all_18S_phy_Vibrio, file = paste0(output_dir, "/Phyloseq_Vi_18S.rds"))