# November 30, 2022
# Same script, just using re-seq files! 

# Hannah - code review comments:
# Hi Joy, this is looking pretty good overall, but there are a some major things I think you should add before you
# move on to downstream analyses. 1) You should run a comparison on an oridnation and with permanova of your blanks and your samples
# to better understand how they're related to each other and if there's any cause to remove your samples if they look like blanks.
# 2) Similarly, to your blanks and samples you should compare your 1st sequencing run to your second sequencing run to make sure that 
# you understand the amount of variation you've introduced by running your samples on separate runs. You can use an ordination and
# permanovas to do this as well. 3) You should be more descriptivie of the "why" in your comments. That way when you revisit this in a couple of months when 
# you're addressing reviewer comments you'll be able to tell them why you made certain decisions. Other than that this is looking pretty 
# good. I've marked all my comments with an "# H:" so you should be able to find them that way.

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
library("ggplot2")      # graphics # H: love that you put what each library is for
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("vegan")        #Used for rarifying

#library(here) # H:adding this in for myself to facilitate loading files, you can comment out, Joy. Similarly , I changed the file paths below.
#here::i_am("reseq_phyloseq_finalobject_script_edit.R")# H:adding this in for myself to facilitate loading files

# Read the data into R from excel file incubation_16S_reseq.xlsx 
seq_tab.raw <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "seqtab_final_reseq") # fill this in with the file path for your DADA2 seqtab output
taxonomy <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "tax_final_reseq", na = c("","NA")) # fill this in with the tax_final file DADA2 output
metadata <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "metadata_final_reseq") # fill this in with the metadata file 


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

# Check to see if taxa were removed
dorm1 # 22923 taxa

# The number of taxa that were removed from being unassigned at the Phylum level: 1249
# Filtering chloroplasts and mitochondria
dorm1 <- subset_taxa(dorm1, Family != "Mitochondria")
dorm1 <- subset_taxa(dorm1, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm1  #13620 taxa
# Number of taxa removed from mitochondria and chloroplasts: 9303 #H: I always like to calculate what this number is as a % of total reads. That way you get a sense of how much non-prokaryotic stuff was in your DNA

# Removing sample blanks from the phyloseq object
dorm_blanks <- subset_samples(dorm1, sample_blank == "BLANK") #making a data frame with just the blanks

# Verify if the blanks contain organisms that are cause for concern/look at the composition of blanks and your actual samples
plot_bar(dorm_blanks, fill = "Phylum") # just the blanks
plot_bar(dorm1, fill = "Phylum") # whole object with blanks included
plot_bar(dorm_noblank, fill = "Phylum") # whole object with no blanks

# Run an nmds of all your blank samples
total_otu_table_blanks <- otu_table(dorm_blanks)

# Transform the OTU table 
dorm_blanks_transformed <- t(sqrt(total_otu_table_blanks))
dorm_blanks_dist <- vegdist(dorm_blanks_transformed, method = "bray")
dorm_blanks.nmds <- metaMDS(dorm_blanks_dist, k = 2, trymax = 1500)

print(dorm_blanks.nmds)
plot(dorm_blanks.nmds)

#*********************************************************************************
# Let's do the same as above for dorm1 (seq, re-seq, and blanks)
total_otu_table_dorm1 <- otu_table(dorm1)
# Transform table
dorm_total_transformed <- t(sqrt(total_otu_table_dorm1))

# Run and NMDS
dorm_total_dist <- vegdist(dorm_total_transformed, method = "bray")
dorm_total.nmds <- metaMDS(dorm_total_transformed, k = 2, trymax = 1500)
# Plot NMDS
print(dorm_total.nmds)
plot(dorm_total.nmds)
#*********************************************************************************
# Now let's do the same for the seq and re-seq data 
# Doing the same as above but with the re-sequenced samples
dorm_reseq <- subset_samples(dorm1, original_reseq == "reseq")# yay this worked! 
reseq_otu_table <- otu_table(dorm_reseq) # need to change this to data frame if you want to do anything with it

# Transform the OTU table
dorm_reseq_transform <- t(sqrt(reseq_otu_table)) 

# Run an NMDS 
dorm_reseq_dist <- vegdist(dorm_reseq_transform, method = "bray") # looking at how similar the re-sequenced samples are
dorm_reseq.nmds <- metaMDS(dorm_reseq_dist, k = 2, trymax = 1500)

plot((dorm_reseq.nmds)) # this makes sense because the samples came from three sites, and they cluster in three groups here
# However there is one extreme outlier
#****************************************************************************************
# Preparing the original sequence data # START HERE!!! JAN 25
orig_seq_counterparts <- read_excel("~/Desktop/incubation_16S_reseq.xlsx", sheet = "orig_seq")
# Somehow create as an OTU table
# Using the following:
library(dplyr)
colnames(df %>% select(where(is.numeric)) %>% select(which(colSums(.) == 0)))

# Transform the OTU table
orig_seq_cp_otu <- as.data.frame(orig_seq_counterparts)
orig_seq_cp <- t(sqrt(orig_seq_cp_otu)) 

otutable <- as.matrix(dorm1@otu_table)
original_seq <- subset(dorm1@otu_table, col()
                       
# Run an NMDS 
orig_seq_cp_dist <- vegdist(dorm_reseq, method = "bray")
dorm_reseq.nmds <- metaMDS(dorm_reseq_dist, k = 2, trymax = 1500)

# Plot NMDS

#Comparing seq originals to re-seq (repeat above)

# Run NMDS and permanova of seq and reseq 
# NMDS
                       
# PERMANOVA 
perm <- adonis2(dorm_reseq_dist ~ seq_dist, permutations = 999, method = "bray")
                       
print(perm)
summary(perm)

#****************************************************************January 25, 2023

# Everything but the reseq samples
seq <- subset.matrix(total_otu_table_dorm1, col = "reseq")
seq.matrix <- as.data.frame(seq) # this contains everything BUT the resequenced data

seq_transformed <- t(sqrt(seq))
seq_dist <- vegdist(seq_transformed, method = "bray")
reseq.nmds <- metaMDS(reseq_dist, k = 2, trymax = 1500)
orditorp(reseq.nmds,display = "sites",cex = 1.25,air = 0.01)
plot(reseq.nmds)

# PERMANOVA of blanks and samples here
blanks.matrix <- as.matrix(total_otu_table_blanks)
total.matrix <- as.matrix(total_otu_table_dorm1)

blanks_samples_aov<- aov(total.matrix ~ blanks.matrix)
print(blanks_samples_aov)
summary(blanks_samples_aov) # I am not sure that this confirms anything--we knew from the beginning that the blanks were different


# Run a PERMANOVA of the seq and re-seq samples here
seq_reseq_aov <- aov(total.matrix ~ reseq)
summary(seq_reseq_aov)
print(seq_reseq_aov)

# Subset the re-sequenced samples and their corresponding originals 
dorm_reseq <- subset(dorm1, sample_ID == "reseq")

dorm_reseqseq <- subset_samples(dorm1, sample == "1, 15, 54, 13, 40, 41, 42, 48, 24, 3, 6, 32, 7, 8, 9, 10")
# Verify if the re-seq samples look okay in terms of comp
plot_bar(dorm_reseq) +
text(dorm_reseq)

# Run an nmds of the samples (reseq and their original)

# Run a Permanova between the re-seq and the corresponding originals?

# Based on the NMDS it looks like the re-seq and the total object are within the same area, which would be indicate little to know variation?
# ****** END COMPARISON OF SEQ, RESEQ, AND BLANKS ***********************************

# Removing the blanks from the samples
dorm_noblank <- subset_samples(dorm1, sample_blank == "Sample")

# Check to see if samples were removed (should be 5 for the re-seq round)
dorm_noblank # this checks out, 5 samples were removed

# Rarifaction step

# Rarify OTU table of the new object (dorm_noblank) #this includes samples that are under 5k 
otu_rarify <- as.data.frame(dorm_noblank@otu_table)
rarecurve(otu_rarify, step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) #H: this is nice! I'll include some code below that I like to use in my rarefaction explorations
#*******************************************************
# H: Hannah's example rarefaction histogram code
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

# Let's rarify at 5k # H: missing in the comments, why did you decide to rarify at 5K?
dorm_noblank.rarefied <- rarefy_even_depth(dorm_noblank, rngseed = 1, sample.size = 5000, replace = FALSE)

# Now let's check the rarifaction plot after rarifying 
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
# H: I agree with Sean! But instead of looking at the samples with <5000 reads in isolation from all the other samples,
# instead, maybe look at them all together (with the blanks) at a very low rarefaction depth. ~1000 or no rarefaction at all.
# The goal of this being to understand how similar your blanks and low read samples are to your high read samples, which you
# can't do if you only look at them in isolation.

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

# End of over5k check #H: missing from your comments here, what did you conclude/find after comparing the under 5k to over 5k rarefaction curves?
#**************************************************************************************

# Make a map data object # I do not recall why this is necessary H: It's just nice to have if you want to plot the NMDS, you have to merge it with the mapping data. So you're extracting the mapping data so you can merge it with the NMDS points and plot them below.
mapdata_noblank <- sample_data(dorm_noblank.rarefied)
mapdata_noblank$rowname <- rownames(mapdata_noblank)
dorm.nmds.mapdata_noblank <- mapdata_noblank$points %>% data.frame() %>% # H: I think there's an error here, mapdata_noblank$points appears to be the output from the function metaMDS(), but you also named your mapping file the same thing. Probably instead of mapdata_noblank, you want to use dorm_noblank.nmds instead; I've experimentally changed it here 
  rownames_to_column(var = "rowname") %>%  # Adds column called "rowname" to metadata and NMDS
  left_join(mapdata_noblank, by = "rowname") %>% # Matches row name from data to sample meta and joins them # H: I think there's a bug in your code here, I think you probably don't need the final "%>%" sign in this part.

# H: Ahh, here it is, nevermind this seems the correct set of code. I recommend deleting lines 175-177 above, to fix it.  
dorm.nmds.mapdata_noblank <- dorm_noblank.nmds$points %>% data.frame() %>% # H: Fixed the indentation here, it was mixed up because of automatic indentation following your dangling pipe (%>%) above.
  rownames_to_column(var = "rowname") %>%   ##adds column called rowname to metadata and NMDS
  left_join(mapdata_noblank, by = "rowname") ## matches row name from data to sample meta and joins them

# Visualize with ggplot2 # H: missing here - commentary on what you're tring to visualize with ggplot. What do you want this particular ordination to show?
ggplot(dorm.nmds.mapdata_noblank, aes(x = MDS1, y = MDS2, color = site)) + #sets the default for the sub statement # H: error, need to remove = sign
  geom_point(aes(shape = pre_post_thaw),alpha = 1, size = 3) +
  geom_line(aes(group = core_rep, color = site))+
  geom_text(aes(label = sample_ID), size = 2, color = "black")

# Note: alpha in aes shows more opaque (seethroughness)
#geom_text(aes(label = sample), size =2)
#***************************************************************************
# END OF SCRIPT 