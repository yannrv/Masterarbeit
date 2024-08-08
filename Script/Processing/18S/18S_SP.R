## Script run on Apphub account of uni Greifswald
## Processing of 'Seastore_project'

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

# Merge multiple runs 18S sequence tables
seqtab_SS1 <- readRDS("/home/rstudio/Downloads/seqtab_18S_ag_t4.rds")
seqtab_SS2 <- readRDS("/home/rstudio/Downloads/seqtab_18S_ag_t5.rds")
seqtab_SS3 <- readRDS("/home/rstudio/Downloads/seqtab_18S_ag.rds")
st.all <- mergeSequenceTables(seqtab_SS1, seqtab_SS2, seqtab_SS3)

# Display the dimensions of the merged sequence table
dim(st.all)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(st.all)))

# Our expected amplicon size is 250 bp, so we discard all samples that are too short or too long, allowing for a 10 bp wiggle room either way
# The amplified V4 region of the 18S rRNA gene does not have great length variability
seqtab_cut <- st.all[,nchar(colnames(st.all)) %in% 240:260]

# Display the dimensions of the trimmed sequence table
dim(seqtab_cut)

# Verify the distribution of sequence lengths after trimming
table(nchar(getSequences(seqtab_cut)))

# Collapse no mismatch sequences with a minimum overlap of 200 bp
seqtab_cut_uni <- collapseNoMismatch(seqtab_cut, minOverlap = 200, verbose = TRUE)

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab_cut_uni, method="consensus", multithread=TRUE, verbose=TRUE)

# Calculate the proportion of chimeras from all reads
chimera_proportion <- sum(seqtab.nochim) / sum(seqtab_cut)
print(chimera_proportion)

# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))  
names(dna) <- DNAStringSet(dna)

# Create a phyloseq object from the sequence (ASV)-abundance-table and the DNA stringset
phy_18S_SP <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), refseq(dna))

# Change ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(phy_18S_SP) <- paste0("ASV", seq(ntaxa(phy_18S_SP)))

# Export ASV and sequence table for classification with CREST4
writeXStringSet(refseq(phy_18S_SP), "Seqs_18S_all.fasta", format = "fasta")
asv_table_df_t_SS <- t(as.data.frame(otu_table(phy_18S_SP)))
write.csv(asv_table_df_t_SS, "Asv_18S_counts_all.csv")

# Save the workspace
save.image("phy_18S_SP.RData")

## After running Crest4 on Jupyter for classification

# Define taxonomic ranks
tax_ranks <- c("ASV", "root", "type", "domain", "rank_I", "rank_II", "rank_III", "rank_IV", "rank_V", "rank_VI", "rank_VII", "rank_VIII", "rank_IX", "rank_X", "rank_XI")

# Read and clean the classification results
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/all_seqs_18S_SS.fasta.crest4/assignments.txt")), 
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
meta <- read.csv("~/nextcloud/New_Meta_SS.csv", header = TRUE, sep = ",", na.strings = c(""," ","NA"), nrows = 128, as.is = TRUE) 
row.names(meta) <- meta$Sample_name_18S

# Replace specific non-ASCII characters in the Site column
meta$Site <- gsub("\xf6", "ö", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe4", "ä", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe5", "å", meta$Site, useBytes = TRUE)

# Merge phyloseq objects
all_18S_phy <- merge_phyloseq(phy_18S_SP, tax_phy, sample_data(meta))
sample_names(all_18S_phy) <- sample_data(all_18S_phy)$Sample_ID

# Print the merged phyloseq object
print(all_18S_phy)

# Save the workspace
save.image("~/Documents/workspace/18S_SS_all.RData")

# Define the output directory
output_dir <- "/home/rstudio/Documents/phyloseq"

# Save the phyloseq object in RDS format
saveRDS(all_18S_phy, file = paste0(output_dir, "/Phyloseq_SS_18S.rds"))