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
dorm <- phyloseq(ASV, TAX, samples) 

# Prune samples from the phyloseq object dorm created above, and call the pruned object dorm1
dorm1 <- prune_samples(sample_sums(dorm) > 0,dorm)
dorm1 #24172 taxa
# Remove taxa that are unassigned at the Phylum level
dorm1 <- subset_taxa(dorm1, !is.na(Phylum)) # Note that sample blanks are still included here 
# The number of taxa that were removed from being unassigned at the Phylum level: 1249
# Filtering chloroplasts and mitochondria
dorm1 <- subset_taxa(dorm1, Family != "Mitochondria")
dorm1 <- subset_taxa(dorm1, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm1  #13620 taxa
# Number of taxa removed from mitochondria and chloroplasts: 9303 # calculate as % reads for publication 

# Removing the blanks from the samples
dorm_noblank <- subset_samples(dorm1, sample_blank == "Sample")

# Check to see if samples were removed (should be 5 for the re-seq round)
dorm_noblank # this checks out, 5 samples were removed

# Rarifaction step

# Joy's rarefaction step

# Rarify OTU table of the new object (dorm_noblank) (I chose to rarefy at 5k because that is where all of the samples are preserved)
otu_rarify <- as.data.frame(dorm_noblank@otu_table)
rarecurve(otu_rarify, step = 50, cex = 0.5, label = FALSE)
dorm_noblank.rarefied.reseq <- rarefy_even_depth(dorm_noblank, rngseed = 1, sample.size = 5000, replace = FALSE)

# Now let's check the rarifaction plot after rarifying 
rarify.check <- as.data.frame(dorm_noblank.rarefied.reseq@otu_table)
rarecurve(rarify.check, step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) 

# Create the rarified OTU table
total_otu_table_noblank_reseq <- otu_table(dorm_noblank.rarefied.reseq)

# Transform the OTU table 
dorm_noblank_transformed_reseq <- t(sqrt(total_otu_table_noblank_reseq))

# Create a distance matrix using vegan and bray
dm_noblank_reseq <- vegdist(dorm_noblank_transformed_reseq, method = "bray")
# I did the following to check how many samples made it through rarifaction
dm_noblank_mat_reseq <- as.matrix(dm_noblank_reseq)

# Run an NMDS (stress values indicate how well the variation is represented, stress less than 0.05 is good, below 0.3 is poor, below 0.2 is okay) pp: NIce comment, here's the link for that info: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
dorm_noblank.nmds <- metaMDS(dm_noblank_reseq, k = 2, trymax = 1500) ## "k" is the number of axes, sometimes 2 is not enough, "trymax" is the number of times it tries to find solutions
print(dorm_noblank.nmds)
plot(dorm_noblank.nmds)
# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) #H: this is nice! I'll in
# Check to see if taxa were removed
dorm1 # 22923 taxa


# Save the phyloseq object to your desktop/folder if this is your final object (save individual parts of the phyloseq object)
saveRDS(dorm_noblank.rarefied.reseq,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_reseq.rds")
saveRDS(dorm_noblank.rarefied.reseq@otu_table,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_OTU_reseq.rds")
saveRDS(dorm_noblank.rarefied.reseq@tax_table,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_tax_table_reseq.rds" )
saveRDS(dorm_noblank.rarefied.reseq@sam_data,"~/GitHub/Masters_publication/Masters_publication_R/dorm1rarefied_sam_data_reseq.rds")

#*********Comparing sequence runs************************************************
# Make a phyloseq object that contains the original and re-sequenced samples
dorm1
# Ensure that the mapping file contains a column that indicates whether or not the sample has been resequenced
# Do all the pre-cleaning steps (normalizing and rarifying, transform OTU table)

perm <- adonis2(dm_noblank_mat_reseq ~ dorm_noblank.rarefied.reseq@sam_data[["original_reseq_neither"]] + dorm_noblank.rarefied.reseq@sam_data[["site"]], permutations = 999)
summary(perm)
print(perm)

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
readcount_hist <- data.frame(ReadCounts = sort(colSums(dorm_noblank@otu_table))) %>%
  rownames_to_column(var = "SampleID") %>%
  # Maybe you want to know if resequenced samples are different read counts; you can either add a column manually like this
  mutate(reseq = ifelse(grepl("reseq", SampleID), "resequenced", "original")) %>%
  # Or you can simply merge your entire mapping data to this data and then the world is your oyster
  left_join(dorm_noblank@sam_data %>% data.frame() %>% rownames_to_column(var = "rownames"), 
            by = c("SampleID" = "rownames")) %>%
  ggplot() +
  geom_histogram(aes(x = ReadCounts)) + 
  #facet_wrap(~reseq) + # H: for example now you can look at the effect of resequencing on read counts;
  facet_wrap(~site) + # H: or see if there are sequencing depth differences between sites;
  geom_vline(xintercept = rar_level, color = "red") +
  ggtitle("Histogram of Read Counts in non-blank samples")
readcount_hist # H: Seems like that highly sequenced tail is all due to your resequnced samples. You'll want to be careful. You may have different communities simply due to the better sequencing depth

# END OF SCRIPT