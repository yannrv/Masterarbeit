# After inferring amplicon sequence variants from the dataset
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

# Merging sequence tables for further processing
seqtab_SS1 <- readRDS("/home/rstudio/Downloads/seqtab_ag_t4.rds")
seqtab_SS2 <- readRDS("/home/rstudio/Downloads/seqtab_16S_ag_t5.rds")
seqtab_SS3 <- readRDS("/home/rstudio/Downloads/seqtab_ag.rds")
seqtab_SS4 <- readRDS("/home/rstudio/Downloads/seqtab_16S_bg_t5.rds")
seqtab_SS5 <- readRDS("/home/rstudio/Downloads/seqtab_bg.rds")
seqtab_SS6 <- readRDS("/home/rstudio/Downloads/seqtab_bg_t4.rds")
st.all <- mergeSequenceTables(seqtab_SS1, seqtab_SS2, seqtab_SS3, seqtab_SS4, seqtab_SS5, seqtab_SS6)

# Checking the dimensions of the merged sequence table
dim(st.all)

# Inspecting the distribution of sequence lengths
table(nchar(getSequences(st.all)))

# Filtering out sequences that are too short or too long, allowing for a 10 bp wiggle room
# The expected amplicon size is 250 bp
seqtab_cut_all <- st.all[, nchar(colnames(st.all)) %in% 240:260]
dim(seqtab_cut_all)

# Inspecting the distribution of sequence lengths after filtering
table(nchar(getSequences(seqtab_cut_all)))

# Collapsing sequences with no mismatches
seqtab_cut_uni_all <- collapseNoMismatch(seqtab_cut_all, minOverlap = 200, verbose = TRUE)

# Removing chimeras
seqtab.nochim_all <- removeBimeraDenovo(seqtab_cut_uni_all, method = "consensus", multithread = TRUE, verbose = TRUE)

# Calculating the proportion of chimeras from all reads
sum(seqtab.nochim_all) / sum(seqtab_cut_all)

# Creating a DNAStringSet from the ASVs
dna_all <- DNAStringSet(getSequences(seqtab.nochim_all))
names(dna_all) <- DNAStringSet(dna_all)

# Creating a phyloseq object from the sequence (ASV)-abundance table and the DNA stringset
all_16S_phy_all <- phyloseq(otu_table(seqtab.nochim_all, taxa_are_rows = FALSE), refseq(dna_all))

# Changing ASV names from the sequence to ASV0001 ... ASV31182
taxa_names(all_16S_phy_all) <- paste0("ASV", seq(ntaxa(all_16S_phy_all)))

# Exporting ASV and sequence table for classification with crest4
writeXStringSet(refseq(all_16S_phy_all), "all_seqs_16S_all.fasta", format = "fasta")
asv_table_df_t_all <- t(as.data.frame(otu_table(all_16S_phy_all)))
write.csv(asv_table_df_t_all, "all_asv_16S_counts_all.csv")

# Saving the R environment
save.image("all_16S_all.RData")

# After running in Python crest4

tax_ranks <- c("ASV", "root", "type", "kingdom", "superkingdom", "superphylum", "phylum", "class", "order", "family", "genus", "species")

# Reading the assignment file and filling NA values
assign <- read.table(text = gsub("\t|; ", ";", readLines("~/all_seqs_16S_all.fasta.crest4/assignments.txt")),
                     header = FALSE, fill = TRUE, na.strings = "", sep = ";", col.names = tax_ranks)

# Filling columns with NA values with the entry from the previous cell and adding "unclassified"
for (i in 1:ncol(assign)) {
  assign[, i] <- ifelse(is.na(assign[, i]), paste0(assign[, i - 1], "_unclassified"), (assign[, i]))
}

# Setting row names
row.names(assign) <- assign$ASV

# Creating a taxonomic matrix
assign_ma <- as.matrix(assign[, -1])
assign_ma <- gsub("(_unclassified_.{1,})", "_unclassified", assign_ma)

# Creating a tax_table object
tax_phy <- tax_table(assign_ma)

# Reading in metadata
meta <- read.csv("~/nextcloud/New_Meta_SS.csv", header = TRUE, sep = ",", na.strings = c("", " ", "NA"), nrows = 128, as.is = TRUE)
row.names(meta) <- meta$Sample_name_16S
meta$Location <- gsub("\xf6", "ö", meta$Location, useBytes = TRUE)
meta$Location <- gsub("\xe4", "ä", meta$Location, useBytes = TRUE)
meta$Location <- gsub("\xe5", "å", meta$Location, useBytes = TRUE)

# Merging phyloseq objects and setting sample names
physeq_16S_SP <- merge_phyloseq(tax_phy, all_16S_phy_all, sample_data(meta))
sample_names(physeq_16S_SP) <- sample_data(all_16S_phy_all)$Sample_ID

# Printing the phyloseq object
print(physeq_16S_SP)

# Saving the R environment
save.image("~/Documents/workspace/16S_SS_all.RData")

# Saving the phyloseq object
output_dir <- "/home/rstudio/Documents/phyloseq"
saveRDS(physeq_16S_SP, file = paste0(output_dir, "/Phyloseq_16S_SP.rds"))
