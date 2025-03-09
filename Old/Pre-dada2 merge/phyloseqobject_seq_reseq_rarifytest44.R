# Joy O'Brien Master's Publication
# March 27, 2023

# Script: This is the beginning of data analysis (resequenced data for publication)
# This script covers the following: 
# load in raw data (sequence table, taxonomy table and sample data from the DADA2 Pipeline), trim data (remove blank controls/mitochondria etc.), rarify data, compare 
# first sequence run to the second round of sequencing
# and create the final phyloseq object that will be used for data analysis
# Made from the reseq_phyloseq_object_FINAL.R script (JMO_masters_2022 repo)



# Install the packages if needed
#install.packages("dplyr")
#install.packages("readxl")    
#install.packages("ggplot2")   
#install.packages("phyloseq")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

# Load packages 
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("vegan")        #Used for rarifying
library("stats")

# Read the data into R from excel file incubation_16S_reseq.xlsx 
seq_tab.raw <- read_excel("~/GitHub/Masters_publication/incubation_16S_reseq.xlsx", sheet = "seqtab_final_reseq") # fill this in with the file path for your DADA2 seqtab output
taxonomy <- read_excel("~/GitHub/Masters_publication/incubation_16S_reseq.xlsx", sheet = "tax_final_reseq", na = c("","NA")) # fill this in with the tax_final file DADA2 output
metadata <- read_excel("~/GitHub/Masters_publication/incubation_16S_reseq.xlsx", sheet = "metadata_final_reseq") # fill this in with the metadata file 

# Define row names from the ASV column 
seq_tab <- seq_tab.raw %>%
  tibble::column_to_rownames("ASV") # H: probably don't need to use tibble:: to call column_to_rownames but it won't hurt.

taxonomy <- taxonomy %>%
  tibble::column_to_rownames("ASV_ID") 

metadata <- metadata %>% 
  tibble::column_to_rownames("sample_name") %>% 
  mutate(sample_blank = ifelse(test = grepl("BLANK", sample_ID), no = "Sample", yes = "BLANK")) %>%
  mutate(core_rep = gsub("_POST_SOIL", "", sample_ID)) %>%
  mutate(core_rep = gsub("_PRE_SOIL", "", core_rep))

# Transform into matrices 
seq_tab <- as.matrix(seq_tab)
taxonomy <- as.matrix(taxonomy)

# Transform to phyloseq objects
ASV <- otu_table(seq_tab, taxa_are_rows = TRUE) 
TAX <- tax_table(taxonomy)
samples = sample_data(metadata)

# Create the phyloseq object
dorm_rt_44 <- phyloseq(ASV, TAX, samples) #dorm_rt = dorm_rarefy test

# Prune samples from the phyloseq object dorm created above, and call the pruned object dorm1
dorm_rt_44 <- prune_samples(sample_sums(dorm_rt_44) > 0,dorm_rt_44)
dorm_rt_44 #24172 taxa
# Remove taxa that are unassigned at the Phylum level
dorm_rt_44 <- subset_taxa(dorm_rt_44, !is.na(Phylum)) # Note that sample blanks are still included here 
# The number of taxa that were removed from being unassigned at the Phylum level: 1249
# Filtering chloroplasts and mitochondria
dorm_rt_44 <- subset_taxa(dorm_rt_44, Family != "Mitochondria")
dorm_rt_44 <- subset_taxa(dorm_rt_44, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm_rt_44 #13620 taxa
# Number of taxa removed from mitochondria and chloroplasts: 9303 # calculate as % reads for publication 

# Removing the blanks from the samples
dorm_rt_44noblank <- subset_samples(dorm_rt_44, sample_blank == "Sample")

# Check to see if samples were removed (should be 5 for the re-seq round)
dorm_rt_44noblank # this checks out, 5 samples were removed

# Rarifaction step

# Joy's rarefaction step

# Rarify OTU table of the new object (dorm_noblank) (I chose to rarefy at 5k because that is where all of the samples are preserved)
otu_rarify <- as.data.frame(dorm_rt_44noblank@otu_table)
rarecurve(otu_rarify, step = 50, cex = 0.5, label = FALSE)
dorm_rt_44noblank.rarefied.reseq <- rarefy_even_depth(dorm_rt_44noblank, rngseed = 1, sample.size = 4400, replace = FALSE)

# Now let's check the rarifaction plot after rarifying 
rarify.check <- as.data.frame(dorm_rt_44noblank.rarefied.reseq@otu_table)
rarecurve(rarify.check, step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) 

# Create the rarified OTU table
total_otu_table_44rt_noblank_reseq <- otu_table(dorm_rt_44noblank.rarefied.reseq)

# Transform the OTU table 
dorm_rt_44noblank_transformed_reseq <- t(sqrt(total_otu_table_44rt_noblank_reseq))

# Create a distance matrix using vegan and bray
dm_rt_44noblank_reseq <- vegdist(dorm_rt_44noblank_transformed_reseq, method = "bray")
# I did the following to check how many samples made it through rarifaction
dm_rt_44noblank_mat_reseq <- as.matrix(dm_rt_44noblank_reseq)

# Run an NMDS (stress values indicate how well the variation is represented, stress less than 0.05 is good, below 0.3 is poor, below 0.2 is okay) pp: NIce comment, here's the link for that info: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
dorm_rt_44noblank.nmds <- metaMDS(dm_rt_44noblank_reseq, k = 2, trymax = 1500) ## "k" is the number of axes, sometimes 2 is not enough, "trymax" is the number of times it tries to find solutions
print(dorm_rt_44noblank.nmds)
plot(dorm_rt_44noblank.nmds)

plot_ordination(
  physeq = dorm_rt_44noblank.rarefied.reseq,
  ordination = dorm_rt_44noblank.nmds,
  color = "site",
  shape = "pre_post_thaw",
  #title = "NMDS of Permafrost Bacterial Communities"
) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(size = 2)+
  geom_text(aes(label = sample), size = 4, color = "black") +
  theme_classic() +
  theme(axis.text = element_text(size = 18)) +
  scale_color_discrete("Site") +
  scale_shape_discrete("Treatment")
scale_shape_manual("Treatment", labels = c("Pre-thaw", "Post-thaw"), values = c(2, 1))

perm <- adonis2(dm_rt_44noblank_mat_reseq ~ dorm_rt_44noblank.rarefied.reseq@sam_data[["original_reseq_neither"]] + dorm_rt_44noblank.rarefied.reseq@sam_data[["site"]], permutations = 999)
summary(perm)
print(perm)
# example: 
perm <- adonis2(dm_noblank_mat_reseq ~ dorm_noblank.rarefied.reseq@sam_data[["original_reseq_neither"]] + dorm_noblank.rarefied.reseq@sam_data[["site"]], permutations = 999)
summary(perm)
print(perm)
# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_rt_noblank)))) #H: this is nice! I'll in

# Save the phyloseq object to your desktop/folder if this is your final object (save individual parts of the phyloseq object)
saveRDS(dorm_rt_noblank.rarefied.reseq,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_rt_reseq.rds")
saveRDS(dorm_rt_noblank.rarefied.reseq@otu_table,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_rt_OTU_reseq.rds")
saveRDS(dorm_rt_noblank.rarefied.reseq@tax_table,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_rt_tax_table_reseq.rds" )
saveRDS(dorm_rt_noblank.rarefied.reseq@sam_data,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_rt_sam_data_reseq.rds")

# Identifying the differences in sample content between the rarefy test set and the 5000 rarefied set

# Convert sam_data elements to character vectors
sam_data1 <- as.character(dorm_rt_noblank.rarefied.reseq@sam_data)
sam_data2 <- as.character(dorm1rarefied_reseq@sam_data)

# Find the differences using setdiff
diff_sam_data <- setdiff(sam_data1, sam_data2)

# Print the differences
print(diff_sam_data)
diff_sam_data <-as.matrix(diff_sam_data)
#*********Comparing sequence runs************************************************
# Make a phyloseq object that contains the original and re-sequenced samples
dorm1
# Ensure that the mapping file contains a column that indicates whether or not the sample has been resequenced
# Do all the pre-cleaning steps (normalizing and rarifying, transform OTU table)

# POST HOC
pairwise_result <- pairwise.adonis2(formula = dm_rt_44noblank_mat_reseq ~ group_variable, permutations = 999)

# Adjust for multiple comparisons (e.g., using Bonferroni correction)
adjusted_pvalues <- p.adjust(pairwise_result$p.values, method = "bonferroni")

# View the results
pairwise_result$comparison <- adjusted_pvalues

# PERMANOVA RESULTS
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = dm_noblank_mat ~ dorm_noblank.rarefied@sam_data[["original_reseq_neither"]] + dorm_noblank.rarefied@sam_data[["site"]], permutations = 999)
#Df SumOfSqs      R2       F Pr(>F)    
#dorm_noblank.rarefied@sam_data[["original_reseq_neither"]]  2   0.9909 0.07208  2.5641  0.003 ** 
###### Attempting betadispr for sequence run variation
##print(beta.seq_reseq)

#________________________________________________________________________________________________________________
#** Extra rarefaction code from HHM**
rar_level <- 5000 # H: tweak this number to see where you can adjust the rarefaction to that makes sense
readcount_hist <- data.frame(ReadCounts = sort(colSums(dorm_rt_noblank@otu_table))) %>%
  rownames_to_column(var = "SampleID") %>%
  # Maybe you want to know if resequenced samples are different read counts; you can either add a column manually like this
  mutate(reseq = ifelse(grepl("reseq", SampleID), "resequenced", "original")) %>%
  # Or you can simply merge your entire mapping data to this data and then the world is your oyster
  left_join(dorm_noblank@sam_data %>% data.frame() %>% rownames_to_column(var = "rownames"), 
            by = c("SampleID" = "rownames")) %>%
  ggplot() +
  geom_histogram(aes(x = ReadCounts)) + 
  facet_wrap(~reseq) + # H: for example now you can look at the effect of resequencing on read counts;
  facet_wrap(~site) + # H: or see if there are sequencing depth differences between sites;
  geom_vline(xintercept = rar_level, color = "red") +
  ggtitle("Histogram of Read Counts in non-blank samples")
readcount_hist # H: Seems like that highly sequenced tail is all due to your resequnced samples. You'll want to be careful. You may have different communities simply due to the better sequencing depth

# END OF SCRIPT