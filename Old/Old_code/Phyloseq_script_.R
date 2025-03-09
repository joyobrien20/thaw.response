# Joy O'Brien
# Masters research 2022
# Phyloseq tutorial(https://vaulot.github.io/tutorials/Phyloseq_tutorial.html) & some lines of code taken from Nate's thesis DNA/RNA script
# 16S permafrost incubation data 
# March 2, 2022

# NOTES/WORKFLOW: 
#1 Rarify OTU table from the phyloseq object 
#2 Visualize rarecurve 
#3 decide how to rarify based on the curve 
#4 transform OTU table / run NMDS code (note to self to modify comments above)
#5 Prepare data for plotting (mapdata)
#6 Visualize via ggplot


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
taxonomy <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "tax_final (2)", na = c("","NA")) # fill this in with the tax_final file (DADA2 output but manipulated in Excel first)
metadata <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final") # fill this in with the metadata file 

## Define row names from the ASV column 
seq_tab <- seq_tab %>%
  tibble::column_to_rownames("ASV")

taxonomy <- taxonomy %>%
  tibble::column_to_rownames("#ASV_ID") 

metadata <- metadata %>% 
  tibble::column_to_rownames("sample_name") %>% 
  mutate(sample_blank = ifelse(test = grepl("BLANK", sample_ID), no = "Sample", yes = "BLANK"))

## Transform into matrices
seq_tab <- as.matrix(seq_tab)
taxonomy <- as.matrix(taxonomy)

## Transform to phyloseq objects
ASV = otu_table(seq_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
samples = sample_data(metadata)

# Create the phyloseq object
dorm <- phyloseq(ASV, TAX, samples) 

# Prune samples from the phyloseq object dorm created above, and call the pruned object dorm1
dorm1 <- prune_samples(sample_sums(dorm) > 0,dorm)
dorm1
# Remove taxa that are unassigned at the Kingdom level  
dorm1 <- subset_taxa(dorm1, !is.na(Phylum)) ##blanks are still included here 

# check to see if taxa were removed
dorm1

dorm1 <- subset_taxa(dorm1, !(Phylum %in% c("", "uncharacterized", "NA")))
# check to see if taxa were removed
dorm1

# Filtering chloroplasts and mitochondria
dorm1 <- subset_taxa(dorm1, Family != "Mitochondria")
dorm1 <- subset_taxa(dorm1, Order != "Chloroplast") 
# check to see if mitochondria and chloroplasts were removed
dorm1 


# Code for subtracting the blank OTUs from the sample OTUs

#dorm_decontam <- subset_taxa(dorm1, otu_table(sample= 154_S58_L002, 155_S59_L002, 156_S60_L002,157_S61_L002, 158_S62_L002, 159_S63_L002, 160_S64_L002, 161_S65_L002, 162_S66_L002, 163_S67_L002, 164_S68_L002, 165_S69_L002, 166_S70_L002, 167_S71_L002, 101_S5_L002, 102_S6_L002, 103_S7_L002, 104_S8_L002, 105_S9_L002, 106_S10_L002, 48_S119_L002, 49_S120_L002, 50_S121_L002 ,51_S122_L002, 52_S123_L002, 53_S124_L002)

#dorm_decontam
#dorm_decontam1<- prune_taxa(taxa_sums(sample=="BLANK"))

#dorm_subract <- prune_taxa(dorm1, )

# Visualize data
sample_names(dorm1)
rank_names(dorm)
sample_variables(dorm)
taxa_names(dorm1)

# Rarify using the vegan package
rarecurve(t(otu_table(dorm1)), step = 50, cex = 0.5, label = FALSE)

# Visualize via histogram
hist(sort(colSums(otu_table(dorm1)))) 

# rarefy without replacement
dorm1.rarefied = rarefy_even_depth(dorm1, rngseed = 1, sample.size = 6000, replace = FALSE) ## removed sample sums but may put back

# Visualize the newly rarified curve
rarecurve(t(otu_table(dorm1.rarefied)), step = 50, cex = 0.5, label = FALSE)

# Create the OTU table with the rarified data
total_otu_table <- otu_table(dorm1.rarefied)

# Calculate relative abundance if appropriate
total = median(sample_sums(dorm1.rarefied))
standf = function(x, t=total) round(t * (x / sum(x)))
dorm1 = transform_sample_counts(dorm1, standf)

# Preparing for ordination
dorm_transformed <- t(sqrt(otu_table(dorm1.rarefied)))## takes sq root of OTU table and does the t command, transform command (changes rows and columns)

# (Optional) Preparing for ordination 
##dorm_transformed <- dorm_transformed[rownames(dorm_transformed)!="insert any samples you want to pluck from your data here"

# Run a distance metric before ordination 
dm <- vegdist(dorm_transformed, method = "bray")

# Run the ordinaton
dorm.nmds <- metaMDS(dm, k = 2, trymax = 1500) ## k is the number of axes, sometimes 2 is not enough, trymax is the number of times it tries to find solutions
## stress values: how well the variation is represented (stress less than 0.05 is good, 0.3 is poor)

# Plot the MDS in vegan
plot(dorm.nmds)
orditorp(dorm.nmds,display = "sites",cex = 1.25,air = 0.01)

# Create a map data object
mapdata <- sample_data(dorm1.rarefied)
mapdata$rowname <- rownames(mapdata)
dorm.nmds.mapdata <- dorm.nmds$points %>% data.frame() %>%
  rownames_to_column(var = "rowname") %>%   # Adds column called "rowname" to metadata and NMDS
  left_join(mapdata, by = "rowname") %>% # Matches row name from data to sample meta and joins them
  
select(sample_ID,
         #starts_with("NMDS1"), starts_with("NMDS2"), ## joined NMDS to mapping data by row name and re-ordered columns
         MDS1, MDS2, everything())
  # Making a "blank" column
  #mutate(sample_blank = ifelse(test = grepl("BLANK", sample_ID), no = "Sample", yes = "BLANK")) %>%

# Get stress value
stress.dorm = paste("stress =", round(dorm.nmds$stress, digits = 4)) ## stress = 0.1105

# Plotting the ordination in ggplot 
ggplot(dorm.nmds.mapdata, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = pre_post_thaw, shape = site),alpha = 1, size = 4) #+  # Alpha in "aes" shows more opaque (seethroughness)
  #geom_text(aes(label = sample), size = 2)

# Performing an ANOVA on the dorm object
dorm.pnova <- adonis(formula = dm ~ sample_blank + pre_post_thaw + site, 
                     data = dorm.nmds.mapdata, 
                     perm = 1000)
# Results of the ANOVA
dorm.pnova

## Analyze via bar plot 
plot_bar(dorm1.rarefied, fill = "Phylum") 
plot_bar(dorm1.rarefied, "sample_ID", fill = "Phylum", facet_grid = ~sample_blank)

# Save the phyloseq object 
saveRDS(dorm1.rarefied,"~/Desktop/dorm1rarefied.rds")

#********************************************************
# Unused code: 

# Remove OTUs that do not appear more than 5 times in more than half the samples. *This removed a ton of taxa for me!!* 
# abundant_dorm_taxa <- genefilter_sample(dorm1, filterfun_sample(function(x) x > 5), A=1) ## decide if you want to do this 
# dorm3 <- prune_taxa(abundant_dorm_taxa, dorm1) ## decide if you want to do this 
# UNUSED CODE THUS FAR: 

# Get blank and sample hulls
# get hulls
# spe.chulls <- plyr::ddply(sb.nmds.mapdata, .(Blank), function(df) df[chull(df$MDS1, df$MDS2), ])
# Plot in ggplot

# ## group_by(Blank) %>%
# mutate(NMDS1BlankAvg = mean(MDS1),
#        NMDS2BlankAvg = mean(MDS2)) %>%
# ungroup() %>%
# group_by(SiteNumber) %>%
# mutate(NMDS1SiteAvg = mean(MDS1),
#        NMDS2SiteAvg = mean(MDS2)) %>%
# ungroup() %>%
# *********************************************************

# ***********************************************************************
# IF YOU NEED TO FILTER OUT THE TAXA FROM THE BLANKS OUT FROM THE SAMPLES
# ***********************************************************************
dorm_subtract <- subset_samples(dorm1, sample_blank == "BLANK")
dorm_subtract


# Vectorized approach, Removing the taxa associated with the blanks from the samples
# step 1, convert to 0 and 1
mat <- otu_table(dorm_subtract)
mat <- as.matrix(data.frame(mat))
mat[mat > 0] <- 1
# Step 2: come up with the row sum 
matrow <- rowSums(mat)
# Remove anything from matrow that is 0

matrow_zero <- matrow[matrow == 0]
# Sort(matrow_nozero)

goodOTUs <- names(matrow_zero)
dorm_taxaremoved <- prune_taxa(goodOTUs, dorm1)

hist(matrow_nozero)

# Save the Phyloseq object if this is your final object (but it won't be)

#************************************************************
# IF YOU NEED TO FILTER OUT THE BLANKS START HERE (removing samples, not taking taxa out)
# ***********************************************************

# STEP 1: Make phyloseq object that includes the blanks (dorm1) you already have dorm1
# Find the blanks in the data and remove from phyloseq object 
dorm_noblank <- subset_samples(dorm1, sample_blank == "Sample")
# check to see if samples were removed (should be 23?)
dorm_noblank

# Rarify OTU table of the new object (dorm_noblank)
rarecurve(t(otu_table(dorm_noblank)), step = 50, cex = 0.5, label = FALSE)

# Visualize with a histogram
hist(sort(colSums(otu_table(dorm_noblank)))) 

# Based on line 85, you will now complete the rarifying process
dorm_noblank.rarefied = rarefy_even_depth(dorm_noblank, rngseed = 1, sample.size = 6000, replace = FALSE) ## removed sample sums but may put back

# Visualize/check rarefaction plot
rarecurve(t(otu_table(dorm_noblank.rarefied)), step = 50, cex = 0.5, label = FALSE)

# Create the rarified OTU table
total_otu_table_noblank <- otu_table(dorm_noblank.rarefied)

# Transform the OTU table 
dorm_noblank_transformed <- t(sqrt(total_otu_table_noblank))

# Create a distance matrix using vegan and bray
dm_noblank <- vegdist(dorm_noblank_transformed, method = "bray")

# Run an NMDS (stress values indicate how well the variation is represented, stress less than 0.05 is good, below 0.3 is poor, below 0.2 is okay)
dorm_noblank.nmds <- metaMDS(dm_noblank, k = 2, trymax = 1500) ## "k" is the number of axes, sometimes 2 is not enough, "trymax" is the number of times it tries to find solutions

# Make a map data object
mapdata_noblank <- sample_data(dorm_noblank.rarefied)
mapdata_noblank$rowname <- rownames(mapdata_noblank)


dorm.nmds.mapdata_noblank <- mapdata_noblank$points %>% data.frame() %>%
  rownames_to_column(var = "rowname") %>%   # Adds column called "rowname" to metadata and NMDS
  left_join(mapdata, by = "rowname") %>% # Matches row name from data to sample meta and joins them

dorm.nmds.mapdata_noblank <- dorm_noblank.nmds$points %>% data.frame() %>%
  rownames_to_column(var = "rowname") %>%   ##adds column called rowname to metadata and NMDS
  left_join(mapdata_noblank, by = "rowname") ## matches row name from data to sample meta and joins them

# Visualize with ggplot2 
ggplot(dorm.nmds.mapdata_noblank, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = pre_post_thaw, shape = soil_biomass_endo_DNA),alpha = 1, size = 3)  ## alpha in aes shows more opaque (seethroughness)
  #geom_text(aes(label = sample), size =2)

# Save the Phyloseq object to your desktop/folder if this is your final object
