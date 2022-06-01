# Blanks filtered from data set script
# Data exploring exercise 
# Joy O'Brien, Ernakovich lab, Masters work 2022
# Date created: February/March 2022, last updated: April 10, 2022
# The purpose of this script is to create a phyloseq object from the files exported from the DADA2 pipeline

# Install the packages needed
#install.packages("dplyr")     # To manipulate data frames
#install.packages("readxl")    # To read Excel files into R
#install.packages("ggplot2")   # for high quality graphics

# Load packages 
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("vegan")        #Used for rarifying


## Read the data into R 
seq_tab <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "seqtab_final_") # fill this in with the file path for your DADA2 seqtab output
taxonomy <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "tax_final", na = c("","NA")) # fill this in with the tax_final file (DADA2 output but manipulated in Excel first)
metadata <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final") # fill this in with the metadata file 

## Define row names from the ASV column 
seq_tab <- seq_tab %>%
  tibble::column_to_rownames("ASV")

taxonomy <- taxonomy %>%
  tibble::column_to_rownames("#ASV_ID") 

metadata <- metadata %>% 
  tibble::column_to_rownames("sample_name") %>% 
  mutate(sample_blank = ifelse(test = grepl("BLANK", sample_ID), no = "Sample", yes = "BLANK")) %>%
  mutate(core_rep = gsub("_POST_SOIL", "", sample_ID)) %>%
  mutate(core_rep = gsub("_PRE_SOIL", "", core_rep))
  
# create a column for the core identifier and the rep identifier 
# Transform into matrices
seq_tab <- as.matrix(seq_tab)
taxonomy <- as.matrix(taxonomy)

## Transform to phyloseq objects
ASV <- otu_table(seq_tab, taxa_are_rows = TRUE) 
TAX <- tax_table(taxonomy)
samples = sample_data(metadata)

# Create the phyloseq object
dorm <- phyloseq(ASV, TAX, samples) 

# Prune samples from the phyloseq object dorm created above, and call the pruned object dorm1
dorm1 <- prune_samples(sample_sums(dorm) > 0,dorm)
dorm1 #9073 taxa
# Remove taxa that are unassigned at the Phylum level
dorm1 <- subset_taxa(dorm1, !is.na(Phylum)) ##blanks are still included here 

# check to see if taxa were removed
dorm1 # 8780 taxa 

# The number of taxa that were removed from being unassigned at the Phylum level: 293 taxa
# Filtering chloroplasts and mitochondria
dorm1 <- subset_taxa(dorm1, Family != "Mitochondria")
dorm1 <- subset_taxa(dorm1, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm1  #5725 taxa

# Number of taxa removed from mitochondria and chloroplasts: 3055

# STEP 1: Make phyloseq object that includes the blanks (dorm1) you already have dorm1
# Find the blanks in the data and remove from phyloseq object 
dorm_noblank <- subset_samples(dorm1, sample_blank == "Sample")
# check to see if samples were removed (should be 23?)
dorm_noblank

# Rarify OTU table of the new object (dorm_noblank)
rarecurve(t(otu_table(dorm_noblank)), step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) 

# You will now complete the rarifying process
dorm_noblank.rarefied = rarefy_even_depth(dorm_noblank, rngseed = 1, sample.size = 5000, replace = FALSE) ## removed sample sums but may put back

# Note for Hannah: I went from 51 samples to 31 samples after removing blanks (approx 5 blanks) 
# and rarifying, is that okay to lost that many samples? --> QQ: How did removing blanks remove samples? If it's because you rarefied and some samples were as low in reads as the blanks, generally those samples tend to be samples that need to be resequenced because they failed.

# Visualize/check rarefaction plot
rarecurve(t(otu_table(dorm_noblank.rarefied)), step = 50, cex = 0.5, label = FALSE)

# Create the rarified OTU table
total_otu_table_noblank <- otu_table(dorm_noblank.rarefied)

# Transform the OTU table 
dorm_noblank_transformed <- t(sqrt(total_otu_table_noblank))

# Create a distance matrix using vegan and bray
dm_noblank <- vegdist(dorm_noblank_transformed, method = "bray")

# Run an NMDS (stress values indicate how well the variation is represented, stress less than 0.05 is good, below 0.3 is poor, below 0.2 is okay) pp: NIce comment, here's the link for that info: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
dorm_noblank.nmds <- metaMDS(dm_noblank, k = 2, trymax = 1500) ## "k" is the number of axes, sometimes 2 is not enough, "trymax" is the number of times it tries to find solutions

# Make a map data object
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
  
  ## alpha in aes shows more opaque (seethroughness)
#geom_text(aes(label = sample), size =2)

summary(dorm_noblank.rarefied)
# Save the phyloseq object to your desktop/folder if this is your final object
saveRDS(dorm_noblank.rarefied,"~/Desktop/dorm1rarefied.rds")
saveRDS(dorm_noblank.rarefied@otu_table,"~/Desktop/dorm1rarefied_OTU.rds")
saveRDS(dorm_noblank.rarefied@tax_table,"~/Desktop/dorm1rarefied_tax_table.rds" )
saveRDS(dorm_noblank.rarefied@sam_data,"~/Desktop/dorm1rarefied_sam_data.rds")


# pp: Love how you saved the RDS object at the end of this. One thing to consider doing as well is to save
# each individual piece of the phyloseq object, OTU table, taxonomy, and metadata that way you can easily hand 
# them off to someone else in spreadsheet form if you need to. This is particularly important for cases like
# publishing data with you paper or thesis.
# General comments: looks great. To take this to the next level, consider creating and saving a table
# of summary stats about your data set and the effects of each filtering step. These are the kinds of things 
# you will report in your thesis/paper. Exampels include: Rarefaction depth, # samples removed, # taxa removed
# # total taxa in final dataset, % mitochondria/chloroplast, etc. 