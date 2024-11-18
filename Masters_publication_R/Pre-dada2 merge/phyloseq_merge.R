# Joy O'Brien Masters Publication
# Merging original and re-sequenced samples that made it through rarefaction at 4400
# Going back to the re-sequeced dataset and going to merge those samples, then rarefy at a new depth
# Then create the final phyloseq object for data analysis

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

# Merge samples that made it past 4400 rarefaction (6 samples)
#*********************************************************************
https://stackoverflow.com/questions/74287436/sam-data-gets-transformed-into-nas-and-numerical-values-when-using-merge-samples 


# Code from stack overflow 
require(phyloseq)
require(tidyverse)

# Concatenate unique values in a vector
concat_unique <- function(vec){
  uniq <- unique(as.character(vec))
  return(paste(uniq, collapse = "/"))
}

# Like psmelt, but only uses the otu_table and sample_data
ps_semi_melt <- function(ps){
  otu_table(ps) %>%
    data.frame(taxid = row.names(.)) %>%
    rename_with(function(x){gsub("X", "", x)}) %>%
    pivot_longer(!taxid, names_to = "sample_id", values_to = "abundance") %>%
    left_join(sample_data(ps) %>%
                data.frame(sample_id = row.names(.)),
              by = "sample_id")
}

# Function that summarizes a vector based on its class
summarise_vec <- function(vec){
  if(class(vec) %in% c("numeric", "integer", "logical")){
    return(mean(vec, na.rm = T))
  } else if (class(vec) %in% c("factor", "character")){
    return(concat_unique(vec))
  } else {
    stop("Error: unknown column type")
  }
}
# Converts a summary df to an otu_table
summ_to_otu_tbl <- function(summ){
  summ %>% 
    select(taxid, sample_id, abundance) %>% 
    pivot_wider(names_from = "sample_id", values_from = "abundance") %>%
    column_to_rownames('taxid') %>%
    as.matrix() %>%
    otu_table(, taxa_are_rows = TRUE)
  }

# Converts a summary df to sample_data
summ_to_sample_dat <- function(summ){
  summ %>% 
    select(!c(taxid, abundance)) %>% 
    unique() %>%
    column_to_rownames('sample_id') %>%
    sample_data()
}

# Function that merges phyloseq samples based on the names of one or more grouping factors
# present in sample_data(ps)
merge_ps_samples <- function(ps, grouping){
  
  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_table(ps) <- phyloseq::otu_table(t(otu_table(ps)), taxa_are_rows = T)
  }
  
  # Convert to long format
  ps_long <- ps_semi_melt(ps)
  
  # Summarise all columns
  summ <- ps_long %>%
    group_by(across(all_of(!!grouping))) %>%
    group_by(taxid, .add = T) %>%
    summarise(across(everything(), summarise_vec)) %>%
    ungroup()
  
  # Convert to otu_table and sample_data
  otu_tbl <- summ_to_otu_tbl(summ)
  sample_dat <- summ_to_sample_dat(summ)
  
  # Create new physeq object
  new_ps <- phyloseq(otu_tbl, sample_dat, tax_table(ps))
  return(new_ps)
}


ps <- dorm_rt_44
merged_ps <- merge_ps_samples(ps, grouping = "sample_ID")

#______________________________________________________________________________________________________


# SHOULD PROBABLY TRY MERGING WITH THE DORM_RT_44NOBLANK DATASET

# Specifying samples to merge 
# samples_to_merge <- (dorm_rt_44noblank, sample == "6", "15", "13", "24", "41", "42")
samples_to_merge <- subset_samples(dorm_rt_44noblank, sample == "6", "15", "13", "24", "41", "42")

merged_physeq <- merge_samples(dorm_rt_44noblank, sample == "6", "15", "13", "24", "41", "42")
print(sample_data(merged_physeq))


# Specify the sample IDs to be merged
#samples_to_merge <- c("AT4_R6_PRE_SOIL", "AT4_R6_POST_SOIL", "FL2C3_R7_POST_SOIL", "FL2C3_R5_POST_SOIL", "BEO2B_R8_POST_SOIL", "BEO1_R5_PRE_SOIL", "BEO1_R6_PRE_SOIL")

# Subset the phyloseq object to include only the specified samples
samples_to_merge <- c("6", "15", "13", "24", "41", "42")
subset_physeq <- subset_samples(dorm_rt_44noblank, sample %in% samples_to_merge)

# Display the sample data before merging (optional)
print(sample_data(subset_physeq))

# sample_data needs to be coerced to a numeric class before merging to avoid NA problem 
df <- as.data.frame(lapply(sample_data(dorm_rt_44noblank),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(dorm_rt_44noblank)
sample_data(dorm_rt_44noblank) <- sample_data(df)

# Merge the specified samples
merged_physeq <- merge_samples(subset_physeq, group = "sample")

# Remove the samples from the phyloseq object if they are there already 
newPhyloObject <- subset_samples(dorm_rt_44noblank, sample != "6" & sample != "15" & sample != "13" & sample != "24" & sample != "41" & sample != "42")
# Merge the merged_physeq with the original physeq

# Code above edited by GPT
# Subset the phyloseq object to include only the specified samples
samples_to_merge <- c("6", "15", "13", "24", "41", "42")
subset_physeq <- subset_samples(dorm_rt_44noblank, sample %in% samples_to_merge)

# Display the sample data before merging (optional)
print(sample_data(subset_physeq))

# Handle NAs and keep factors intact
df <- as.data.frame(lapply(sample_data(dorm_rt_44noblank), function(y) if (class(y) != "factor") as.factor(y) else y), stringsAsFactors = TRUE)
row.names(df) <- sample_names(dorm_rt_44noblank)
sample_data(dorm_rt_44noblank) <- sample_data(df)

# Merge the specified samples
merged_physeq <- merge_samples(subset_physeq, group = "sample")

# Remove the samples from the original phyloseq object
newPhyloObject <- subset_samples(dorm_rt_44noblank, !(sample %in% samples_to_merge))

# Merge the merged_physeq with the original physeq
final_physeq <- merge_phyloseq(newPhyloObject, merged_physeq)


# Almost there, they are merging but the sample data isn't merging correctly and it's causing NAs and different 
# sample names in the final object


# Display the sample data after merging (optional)
print(sample_data(final_physeq))
merge_trial <- as.character((sample_data(dorm_rt_44noblank)))
merged_physeq <- transform_sample_data(merged_physeq, as.character)

# seems like it goes in the order of subsetting the samples to be merged, merging them, and then adding them back into the phyloseq object? 
# pick up here Jan. 5 




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