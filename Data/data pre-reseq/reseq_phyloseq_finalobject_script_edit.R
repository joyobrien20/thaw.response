# November 30, 2022
# Same script, just using re-seq files! 

# Joy O'Brien
# Master's project, Ernakovich lab 
# Blanks filtered from data set script (meaning that blanks are REMOVED from the final phyloseq object)
# Date created: February/March 2022, last updated: April 10, 2022 & November 30, 2022 for re-seq version
# This is the final version of the script 

# The purpose of this script is to create a phyloseq object from the files exported from the DADA2 pipeline

# Install the packages if needed
#install.packages("dplyr")     # To manipulate data frames
#install.packages("readxl")    # To read Excel files into R
#install.packages("ggplot2")   # for high quality graphics
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
seq_tab <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "seqtab_final_reseq") # fill this in with the file path for your DADA2 seqtab output
taxonomy <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "tax_final_reseq", na = c("","NA")) # fill this in with the tax_final file DADA2 output
metadata <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "metadata_final_reseq") # fill this in with the metadata file 

# Define row names from the ASV column 
seq_tab <- seq_tab %>%
  tibble::column_to_rownames("ASV")

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

# Check to see if taxa were removed
dorm1 # 22923 taxa

# The number of taxa that were removed from being unassigned at the Phylum level: 1249
# Filtering chloroplasts and mitochondria
dorm1 <- subset_taxa(dorm1, Family != "Mitochondria")
dorm1 <- subset_taxa(dorm1, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm1  #13620 taxa
# Number of taxa removed from mitochondria and chloroplasts: 9303

# Removing sample blanks from the phyloseq object
# Find the blanks in the data and remove from phyloseq object 
dorm_noblank <- subset_samples(dorm1, sample_blank == "Sample")

# Check to see if samples were removed (should be 5 for the re-seq round)
dorm_noblank # this checks out, 5 samples were removed

# Rarifaction step

# Rarify OTU table of the new object (dorm_noblank) #this includes samples that are under 5k 
otu_rarify <- as.data.frame(dorm_noblank@otu_table)
rarecurve(otu_rarify, step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank))))

# Let's rarify at 5k 
dorm_noblank.rarefied <- rarefy_even_depth(dorm_noblank, rngseed = 1, sample.size = 5000, replace = FALSE)

# Now let's check the rarifation plot after rarifying 
rarify.check <- as.data.frame(dorm_noblank.rarefied@otu_table)
rarecurve(rarify.check, step = 50, cex = 0.5, label = FALSE)

# Create the rarified OTU table
total_otu_table_noblank <- otu_table(dorm_noblank.rarefied)

# Transform the OTU table 
dorm_noblank_transformed <- t(sqrt(total_otu_table_noblank))

# Create a distance matrix using vegan and bray
dm_noblank <- vegdist(dorm_noblank_transformed, method = "bray")

# Run an NMDS (stress values indicate how well the variation is represented, stress less than 0.05 is good, below 0.3 is poor, below 0.2 is okay) pp: NIce comment, here's the link for that info: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
dorm_noblank.nmds <- metaMDS(dm_noblank, k = 2, trymax = 1500) ## "k" is the number of axes, sometimes 2 is not enough, "trymax" is the number of times it tries to find solutions
print(dorm_noblank.nmds)


# Save the phyloseq object to your desktop/folder if this is your final object (save individual parts of the phyloseq object)
saveRDS(dorm_noblank.rarefied,"~/Desktop/dorm1rarefied_reseq.rds")
saveRDS(dorm_noblank.rarefied@otu_table,"~/Desktop/dorm1rarefied_OTU_reseq.rds")
saveRDS(dorm_noblank.rarefied@tax_table,"~/Desktop/dorm1rarefied_tax_table_reseq.rds" )
saveRDS(dorm_noblank.rarefied@sam_data,"~/Desktop/dorm1rarefied_sam_data_reseq.rds")


# End phyloseq code 
#**********************************************************************************************************************
# Testing to see whether removing samples under 5000 reads from the phyloseq object prior to rarifaction changes anything 
# Note from Sean regarding rarifaction
# So the generating of the rarefaction plot is kind of like a branching point +
# Once you have that you go back to the phyloseq object to actual ratify

# Pretty much you do import it to phyloseq +
# do all your cleaning then take an otu out of phyloseq
# do your rarefaction plots and use that info to ratify to min seq depth
#code from sean

# Let's look at the samples that are under 5k
under_5k <- prune_samples(sample_sums(dorm_noblank)<=5000, dorm_noblank)

# List of samples have under 5k reads (majority of these samples were re-sequenced! and none of the re-sequenced samples are under 5k, so this is good)
sample_names(under_5k)
sample_sum_df <- data.frame(sum = sample_sums(under_5k))

over_5k <- prune_samples(sample_sums(dorm_noblank) >= 5000,dorm_noblank) # this should include only samples that are larger than 5000

# This is what we compare to the rarifaction curve above (otu_rarify) as that one contains all samples (above and below 5k)
otu_over5k <- as.data.frame(over_5k@otu_table)
rarecurve(otu_over5k, step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(over_5k))))

# Let's rarify with this 
dorm_noblank_over5k.rarified <-  rarefy_even_depth(over_5k, rngseed = 1, sample.size = 5000, replace = FALSE)

# Now let's check the rarifation plot after rarifying 
rarify_over5k.check <- as.data.frame(dorm_noblank_over5k.rarified@otu_table)
rarecurve(rarify_over5k.check, step = 50, cex = 0.5, label = FALSE) # the plot doesnt change as the one above-and I am not sure that it's supposed to

# End of over5k check 
#**************************************************************************************

# Make a map data object # I do not recall why this is necessary 
mapdata_noblank <- sample_data(dorm_noblank.rarefied)
mapdata_noblank$rowname <- rownames(mapdata_noblank)
dorm.nmds.mapdata_noblank <- mapdata_noblank$points %>% data.frame() %>%
  rownames_to_column(var = "rowname") %>%  # Adds column called "rowname" to metadata and NMDS
  left_join(mapdata_noblank, by = "rowname") %>% # Matches row name from data to sample meta and joins them 
  
  dorm.nmds.mapdata_noblank <- dorm_noblank.nmds$points %>% data.frame() %>% 
  rownames_to_column(var = "rowname") %>%   ##adds column called rowname to metadata and NMDS
  left_join(mapdata_noblank, by = "rowname") ## matches row name from data to sample meta and joins them

# Visualize with ggplot2 
ggplot(dorm.nmds.mapdata_noblank=, aes(x = MDS1, y = MDS2, color = site)) + #sets the default for the sub statement
  geom_point(aes(shape = pre_post_thaw),alpha = 1, size = 3) +
  geom_line(aes(group = core_rep, color = site))+
  geom_text(aes(label = sample_ID), size = 2, color = "black")

# Note: alpha in aes shows more opaque (seethroughness)
#geom_text(aes(label = sample), size =2)
#***************************************************************************
# END OF SCRIPT 