#############################
# 18S rRNA Phyloseq Pipeline
#############################

# ===== 1. LOAD LIBRARIES AND SET RANDOM SEED =====

# Load required packages
library(dada2)        
library(DECIPHER)     
library(phyloseq)     
library(speedyseq)    
library(ggplot2)      
library(Biostrings)   
library(ape)          
library(microbiome)   
library(vegan)        

# Set seed for reproducibility
set.seed(68)

# ===== 2. IMPORT AND MERGE SEQUENCE TABLES =====

# Import sequence tables 
seqtab_SS1 <- readRDS("~/Documents/Seqtab/seqtab_18S_ag_t4.rds")
seqtab_SS2 <- readRDS("~/Documents/Seqtab/seqtab_18S_ag_t5.rds")
seqtab_SS3 <- readRDS("~/Documents/Seqtab/seqtab_18S_ag.rds")
seqtab_SS4 <- readRDS("~/Documents/Seqtab/seqtab_18S_bg_t4.rds")
seqtab_SS5 <- readRDS("~/Documents/Seqtab/seqtab_18S_bg_t5.rds")
seqtab_SS6 <- readRDS("~/Documents/Seqtab/seqtab_18S_bg.rds")
seqtab_SS7 <- readRDS("~/Documents/Seqtab/Seqtab_18S_Germ.rds")
seqtab_SS8 <- readRDS("~/Documents/Seqtab/Seqtab_18S_HS.rds")
seqtab_SS9 <- readRDS("~/Documents/Seqtab/Seqtab_18S_Salinity.rds")
seqtab_SS10 <- readRDS("~/Documents/Seqtab/Seqtab_18S_V.rds")

# Merge all sequence tables 
st.all <- mergeSequenceTables(seqtab_SS1, seqtab_SS2, seqtab_SS3, seqtab_SS4, 
                              seqtab_SS5, seqtab_SS6, seqtab_SS7, seqtab_SS8, 
                              seqtab_SS9, seqtab_SS10)

#Check dimension
dim(st.all)

# ===== 3. PROCESS AND FILTER SEQUENCE DATA =====

# Inspect distribution of sequence lengths
table(nchar(getSequences(st.all)))

# Collapse highly similar sequences (200bp minimum overlap)
seqtab_uni <- collapseNoMismatch(st.all, minOverlap = 200, verbose = TRUE)


# Remove chimeric sequences 
seqtab.nochim <- removeBimeraDenovo(seqtab_uni, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)


# Calculate proportion of non-chimeric reads
non_chimeric_proportion <- sum(seqtab.nochim) / sum(seqtab_uni)
print(non_chimeric_proportion)


# ===== 4. CREATE INITIAL PHYLOSEQ OBJECT =====

# Create a DNAStringSet from the ASVs for reference
dna <- DNAStringSet(getSequences(seqtab.nochim))
names(dna) <- DNAStringSet(dna)

# Create phyloseq object from sequence table and reference sequences
all_18S_phy_all <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                            refseq(dna))

# Rename ASVs 
taxa_names(all_18S_phy_all) <- paste0("ASV", seq(ntaxa(all_18S_phy_all)))

# Export ASV sequences for taxonomic classification with CREST4
writeXStringSet(refseq(all_18S_phy_all), "all_seqs_18S_Study.fasta", format = "fasta")

# Export ASV count table
asv_table_df_t_SS <- t(as.data.frame(otu_table(all_18S_phy_all)))
write.csv(asv_table_df_t_SS, "all_asv_18S_counts_Study.csv")

# Save workspace image
save.image("all_18S_Study.RData")

# ===== 5. TAXONOMIC ASSIGNMENT WITH CREST4 =====
## (CREST4 is run externally)

# Define taxonomic ranks for CREST4 output
tax_ranks <- c("ASV", "root", "type", "domain", "rank_I", "rank_II", "rank_III", 
               "rank_IV", "rank_V", "rank_VI", "rank_VII", "rank_VIII",
               "rank_IX", "rank_X", "rank_XI")

# Import and process taxonomy assignments from CREST4
assign <- read.table(
  text = gsub("\t|; ", ";", 
              readLines("~/Documents/fasta_files/all_seqs_18S_Study.fasta.crest4/assignments.txt")), 
  header = FALSE, fill = TRUE, na.strings = "", sep = ";", col.names = tax_ranks
)

# Fill in NA values 
for (i in 1:ncol(assign)) {
  assign[, i] <- ifelse(is.na(assign[, i]), 
                        paste0(assign[, i-1], "_unclassified"), 
                        (assign[, i]))
}

# Prepare taxonomy table
row.names(assign) <- assign$ASV
assign_ma <- as.matrix(assign[, -1])  # Remove ASV column
# Clean up multiple occurrences of "_unclassified"
assign_ma <- gsub("(_unclassified_.{1,})", "_unclassified", assign_ma)

# Create taxonomy table for phyloseq
tax_phy <- tax_table(assign_ma)

# ===== 6. IMPORT METADATA AND CREATE COMPLETE PHYLOSEQ OBJECT =====

# Import metadata
meta <- read.csv("~/Metadata/Metadata.Mix.csv", header = TRUE, sep = ",", 
                 na.strings = c("", " ", "NA"), nrows = 241, as.is = TRUE)

# Fix special characters
row.names(meta) <- meta$Sample_name_18S
meta$Site <- gsub("\xf6", "ö", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe4", "ä", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe5", "å", meta$Site, useBytes = TRUE)

# Merge phyloseq components: ASVs, taxonomy, and metadata
all_18S_phy <- merge_phyloseq(all_18S_phy_all, tax_phy, sample_data(meta))
sample_names(all_18S_phy) <- sample_data(all_18S_phy)$Sample_ID

# Print summary of phyloseq object
print(all_18S_phy)

# Save workspace and phyloseq object
save.image("~/Documents/workspace/18S_SS_mix.RData")
output_dir <- "/home/rstudio/Documents/phyloseq"
saveRDS(all_18S_phy, "/Phyloseq_SS_18S_mix.rds")


# ===== 7. FILTER UNNECESSARY TAXA =====

# Remove 
clean_18S_mg <- subset_taxa(all_18S_phy, 
                            !(root == "No hits" | 
                                type %in% c("Archaea", "root_unclassified") | 
                                rank_III == "Embryophyta"))


# Proportion of reads retained
read_proportion <- sum(otu_table(clean_18S_mg)) / sum(otu_table(all_18S_phy))
print(read_proportion)

# Proportion of ASVs retained
asv_proportion <- ncol(otu_table(clean_18S_mg)) / ncol(otu_table(all_18S_phy))
print(asv_proportion)

# Add library size to sample data
sample_data(clean_18S_mg)$LibrarySize_18S <- sample_sums(clean_18S_mg)

# Check library size
print(max(sample_data(clean_18S_mg)$LibrarySize_18S))
print(mean(sample_data(clean_18S_mg)$LibrarySize_18S))
print(min(sample_data(clean_18S_mg)$LibrarySize_18S))


# plot of library sizes
plot(sample_data(clean_18S_mg)$LibrarySize_18S, 
     main = "18S Library Size Distribution",
     xlab = "Sample Index", ylab = "Read Count")
abline(1000, 0, col = "red", lty = 2)  


# ===== 8. FILTER SAMPLES AND PERFORM RAREFACTION =====

# Remove samples with fewer than 1000 reads
clean_18S_mg.2 <- prune_samples(sample_sums(clean_18S_mg) >= 1000, clean_18S_mg)
cat("Samples retained after minimum read filter:", nsamples(clean_18S_mg.2), 
    "of", nsamples(clean_18S_mg), "\n")

# Perform rarefaction 
clean_18S_mg_rarefied <- rarefy_even_depth(
  clean_18S_mg.2, 
  sample.size = min(sample_sums(clean_18S_mg.2)), 
  rngseed = 100, 
  replace = FALSE
)


# ===== 9. CALCULATE ALPHA DIVERSITY METRICS =====

# Calculate alpha diversity metrics (Observed, Shannon, and Pielou's evenness)
alpha_div <- estimate_richness(clean_18S_mg_rarefied, 
                               measures = c("Observed", "Shannon", "Pielou"))

# Add the metrics to the sample data
sample_data(clean_18S_mg_rarefied)$Richness_18S <- alpha_div$Observed
sample_data(clean_18S_mg_rarefied)$Shannon_18S <- alpha_div$Shannon
sample_data(clean_18S_mg_rarefied)$Evenness_18S <- alpha_div$Pielou

# Check results
print(summary(alpha_div))

# ===== 10. EXPORT FINAL DATA AND SAVE =====

# Extract OTU table, sample data, and taxonomy table
otu_df <- as.data.frame(otu_table(clean_18S_mg_rarefied))
sample_df <- as.data.frame(sample_data(clean_18S_mg_rarefied))
tax_df <- as.data.frame(tax_table(clean_18S_mg_rarefied))

# Save the extracted tables as CSV files
write.csv(otu_df, file = "ASV_tab_18S_mix.csv", row.names = TRUE)
write.csv(tax_df, file = "Taxonomy_tab_18S_mix.csv", row.names = TRUE) 
write.csv(sample_df, file = "Sample_data_18S_mix.csv", row.names = FALSE)

# Save the final cleaned and rarefied phyloseq object
saveRDS(clean_18S_mg_rarefied, file = paste0(output_dir, "/Phyloseq_SS_18S_mix_clean.rds"))
