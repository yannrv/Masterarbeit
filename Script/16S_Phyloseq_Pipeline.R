#############################
# 16S rRNA Phyloseq Pipeline
#############################

# ===== 1. LOAD LIBRARIES AND SET RANDOM SEED =====

# Load required packages for 16S rRNA analysis
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
set.seed(100)

# ===== 2. IMPORT AND MERGE SEQUENCE TABLES =====

# Import sequence tables 
seqtab_SS1 <- readRDS("seqtab_ag_t4.rds")
seqtab_SS2 <- readRDS("seqtab_16S_ag_t5.rds")
seqtab_SS3 <- readRDS("seqtab_ag.rds")
seqtab_SS4 <- readRDS("seqtab_16S_bg_t5.rds")
seqtab_SS5 <- readRDS("seqtab_bg.rds")
seqtab_SS6 <- readRDS("seqtab_bg_t4.rds")
seqtab_HS <- readRDS("Seqtab_16S_HS.rds")
seqtab_Sal <- readRDS("Seqtab_16S_Salinity.rds")
seqtab_V <- readRDS("Seqtab_16S_V.rds")
seqtab_Germ <- readRDS("Seqtab_16S_Germ.rds")

# Merge all sequence tables 
st.all <- mergeSequenceTables(seqtab_SS1, seqtab_SS2, seqtab_SS3, 
                              seqtab_SS4, seqtab_SS5, seqtab_SS6,
                              seqtab_HS, seqtab_Sal, seqtab_V, seqtab_Germ)

#Check dimension
dim(st.all)

# ===== 3. PROCESS AND FILTER SEQUENCE DATA =====

# Filter sequences by length (240-260 bp for V4 region of 16S)
seqtab_cut_all <- st.all[,nchar(colnames(st.all)) %in% 240:260]
print(ncol(seqtab_cut_all))

# Collapse identical sequences with minimum 200bp overlap
seqtab_cut_uni_all <- collapseNoMismatch(seqtab_cut_all, minOverlap = 200, verbose = TRUE)
print(ncol(seqtab_cut_uni_all))

# Remove chimeric sequences 
seqtab.nochim_all <- removeBimeraDenovo(seqtab_cut_uni_all, method = "consensus", 
                                        multithread = TRUE, verbose = TRUE)

# Calculate proportion of non-chimeric reads
chimera_ratio <- sum(seqtab.nochim_all) / sum(seqtab_cut_uni_all)
cat("Proportion of reads retained after chimera removal:", chimera_ratio, "\n")
cat("Number of sequences after chimera removal:", ncol(seqtab.nochim_all), "\n")

# ===== 4. CREATE INITIAL PHYLOSEQ OBJECT =====

# Create a DNAStringSet from the ASVs for reference
dna_all <- DNAStringSet(getSequences(seqtab.nochim_all))
names(dna_all) <- DNAStringSet(dna_all)

# Create phyloseq object from sequence table and reference sequences
all_16S_phy_all <- phyloseq(otu_table(seqtab.nochim_all, taxa_are_rows = FALSE), 
                            refseq(dna_all))

# Rename ASVs 
taxa_names(all_16S_phy_all) <- paste0("ASV", seq(ntaxa(all_16S_phy_all)))

# Export ASV sequences for taxonomic classification with CREST4
writeXStringSet(refseq(all_16S_phy_all), "all_seqs_16S_study.fasta", format = "fasta")

# Export ASV count table
asv_table_df_t_all <- t(as.data.frame(otu_table(all_16S_phy_all)))
write.csv(asv_table_df_t_all, "all_asv_16S_counts_study.csv")

# ===== 5. IMPORT AND PROCESS TAXONOMIC ASSIGNMENTS =====
## (CREST4 is run externally)

# Define taxonomic ranks for CREST4 output
tax_ranks <- c("ASV", "root", "type", "kingdom", "superkingdom", "superphylum", 
               "phylum", "class", "order", "family", "genus", "species")

# Import and process taxonomy assignments from CREST4
assign <- read.table(
  text = gsub("\t|; ", ";", 
              readLines("all_seqs_16S_study.fasta.crest4/assignments.txt")), 
  header = FALSE, fill = TRUE, na.strings = "", sep = ";", 
  col.names = tax_ranks
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
meta <- read.csv("Metadata_mix.csv", header = TRUE, sep = ",", 
                 na.strings = c("", " ", "NA"), nrows = 241, as.is = TRUE)
row.names(meta) <- meta$Sample_name_16S

# Fix special characters 
meta$Site <- gsub("\xf6", "ö", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe4", "ä", meta$Site, useBytes = TRUE)
meta$Site <- gsub("\xe5", "å", meta$Site, useBytes = TRUE)

# Create complete phyloseq object with taxonomy, sequences, and metadata
phy_all <- merge_phyloseq(tax_phy, all_16S_phy_all, sample_data(meta))
sample_names(phy_all) <- sample_data(phy_all)$Sample_ID

# Print summary of complete phyloseq object
print(phy_all)

# ===== 7. FILTER UNNECESSARY TAXA =====

# Remove 
filtered_obj <- subset_taxa(phy_all, 
                            !(type %in% c("Chloroplast", "Mitochondria", 
                                          "root_unclassified") | 
                                root == "No hits"))

# Keep only main bacterial genome sequences
clean_16S_mg <- subset_taxa(filtered_obj, type == "Main genome")

# Proportion of reads retained
read_proportion <- sum(otu_table(clean_16S_mg)) / sum(otu_table(phy_all))
print(read_proportion)

# Proportion of ASVs retained
asv_proportion <- ncol(otu_table(clean_16S_mg)) / ncol(otu_table(phy_all))
print(asv_proportion)


# Add library size to sample data
sample_data(clean_16S_mg)$LibrarySize_16S <- sample_sums(clean_16S_mg)

# Check library size 
print(max(sample_data(clean_16S_mg)$LibrarySize_16S))
print(mean(sample_data(clean_16S_mg)$LibrarySize_16S))
print(min(sample_data(clean_16S_mg)$LibrarySize_16S))


# plot of library sizes
plot(sample_data(clean_16S_mg)$LibrarySize_16S, 
     main = "16S Library Size Distribution",
     xlab = "Sample Index", ylab = "Read Count")
abline(4000, 0, col = "red", lty = 2)  



# ===== 8. PERFORM RAREFACTION =====

# Remove samples with fewer than 4000 reads
clean_16S_mg.2 <- prune_samples(sample_sums(clean_16S_mg) >= 4000, clean_16S_mg)
cat("Samples retained after minimum read filter:", nsamples(clean_16S_mg.2), 
    "of", nsamples(clean_16S_mg), "\n")

# Perform rarefaction 
clean_16S_mg_rarefied <- rarefy_even_depth(
  clean_16S_mg.2, 
  sample.size = min(sample_sums(clean_16S_mg.2)), 
  rngseed = 100, 
  replace = FALSE
)


# ===== 9. CALCULATE ALPHA DIVERSITY METRICS =====

# Calculate alpha diversity metrics (Observed, Shannon, and Pielou's evenness)
alpha_div <- estimate_richness(clean_16S_mg_rarefied, 
                               measures = c("Observed", "Shannon", "Pielou"))

# Add the metrics to the sample data
sample_data(clean_16S_mg_rarefied)$Richness_16S <- alpha_div$Observed
sample_data(clean_16S_mg_rarefied)$Shannon_16S <- alpha_div$Shannon
sample_data(clean_16S_mg_rarefied)$Evenness_16S <- alpha_div$Pielou

# Check results
print(summary(alpha_div))


# ===== 10. EXPORT FINAL DATA AND SAVE =====

# Extract OTU table, sample data, and taxonomy table
otu_table <- as.data.frame(otu_table(clean_16S_mg_rarefied))
sample_data <- as.data.frame(sample_data(clean_16S_mg_rarefied))
taxonomy_table <- as.data.frame(tax_table(clean_16S_mg_rarefied))

# Save the extracted tables as CSV files
write.csv(otu_table, "ASV_tab_16S_mix.csv", row.names = TRUE)
write.csv(taxonomy_table, "Taxonomy_tab_16S_mix.csv", row.names = TRUE)
write.csv(sample_data, "Sample_data_16S_mix.csv", row.names = FALSE)

# Save the final cleaned and rarefied phyloseq object
saveRDS(clean_16S_mg_rarefied, "Phyloseq_16S_mix_clean.rds")
