#### Script run on Apphub account of uni Greifswald
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

# Load the 18S sequence table
seqtab_Germ <- readRDS("/home/rstudio/Documents/Seqtab/Seqtab_18S_Germ.rds")

# Display the dimensions of the sequence table
dim(seqtab_Germ)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab_Germ)))

# Our expected amplicon size is 250 bp, so we discard all samples that are too short or too long, allowing for a 10 bp wiggle room either way
# The amplified V4 region of the 18S rRNA gene does not have great length variability
seqtab_cut_Germ <- seqtab_Germ[,nchar(colnames(seqtab_Germ)) %in% 240:260]

# Display the dimensions of the trimmed sequence table
dim(seqtab_cut_Germ)

# Verify the distribution of sequence lengths after trimming
table(nchar(getSequences(seqtab_cut_Germ)))

# Collapse no mismatch sequences with a minimum overlap of 200 bp
seqtab_cut_uni_Germ <- collapseNoMismatch(seqtab_cut_Germ, minOverlap = 200, verbose = TRUE)

# Remove chimeric sequences
seqtab.nochim_Germ <- removeBimeraDenovo(seqtab_cut_uni_Germ, method="consensus", multithread=TRUE, verbose=TRUE)

# Create a DNAStringSet from the ASVs
dna_Germ <- DNAStringSet(getSequences(seqtab.nochim_Germ))  
names(dna_Germ) <- DNAStringSet(dna_Germ)

# Create a phyloseq object from the sequence (ASV)-abundance-table and the DNA stringset
all_18S_phy_Germ <- phyloseq(otu_table(seqtab.nochim_Germ, taxa_are_rows=FALSE), refseq(dna_Germ))

# Change ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(all_18S_phy_Germ) <- paste0("ASV", seq(ntaxa(all_18S_phy_Germ)))

# Export ASV and sequence table for classification with CREST4
writeXStringSet(refseq(all_18S_phy_Germ), "all_seqs_18S_Germ.fasta", format = "fasta")
asv_table_df_t_Germ <- t(as.data.frame(otu_table(all_18S_phy_Germ)))
write.csv(asv_table_df_t_Germ, "all_18S_counts_Germ.csv")

# Save the workspace
save.image("all_18S_Germ.RData")

## After running CREST4 for classification

# Define taxonomic ranks
tax_ranks <- c("ASV", "root", "type", "domain", "rank_I", "rank_II", "rank_III", "rank_IV", "rank_V", "rank_VI", "rank_VII", "rank_VIII",
               "rank_IX", "rank_X", "rank_XI")

# Read and clean the classification results
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/all_seqs_18S_Germ.fasta.crest4/assignments.txt")), 
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
pro_otu <- read.csv("~/all_18S_counts_Germ.csv", header = TRUE, row.names = 1, sep= ",")

# Transpose the OTU table to put sample names in rows
pro_otu <- t(pro_otu)

# Extract sample names
sample_names <- rownames(pro_otu)

# Remove the prefix "X" from sample names
sample_names <- sub("^X", "", sample_names)

# Replace dots with hyphens in sample names
sample_names <- gsub("\\.", "-", sample_names)

# Extract short sample names
sample_names <- sapply(strsplit(sample_names, "1183F-1443R-P01-"), `[`, 2)
short_names <- as.vector(t(sample_names))
row.names(pro_otu) <- short_names

# Read in metadata
meta <- read.csv("~/Metadata/New_Meta_Germ.csv", header = TRUE, sep = ",", na.strings = c(""," ","NA"), nrows = 9, as.is = TRUE) 
row.names(meta) <- meta$Sample_names

# Create an OTU table
pro_OTU <- otu_table(pro_otu, taxa_are_rows=FALSE)

# Merge phyloseq objects
phys_Germ <- merge_phyloseq(pro_OTU, tax_phy, all_18S_phy_Germ, sample_data(meta))
sample_names(phys_Germ) <- sample_data(phys_Germ)$Sample_names
print(phys_Germ)

# Define the output directory
output_dir <- "/home/rstudio/Documents/phyloseq"

# Save the phyloseq object in RDS format
saveRDS(phys_Germ, file = paste0(output_dir, "/Phyloseq_Germ_18S.rds"))