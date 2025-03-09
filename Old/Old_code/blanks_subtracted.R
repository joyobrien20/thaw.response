# Blank taxa subtracted from sample taxa
# Data exploring exercise
# Joy O'Brien, Ernakovich lab, Masters 2022

# Tell R to subtract the blanks from the samples (in terms of ASV, remove ASV in the blanks from the samples)
# NOTE: THIS IS JUST A DATA EXPLORING EXERCISE AND WAS NOT USED IN THESIS/PUBLICATION

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

## Read the data into R (EXCEL SHEET AVAILABLE IN ERNAKOVICH LAB BOX UNDER OBRIEN THESIS)
seq_tab <- read_excel("~/Desktop/incubation_16S_final.xlsx", sheet = "seqtab_final_") # fill this in with the file path for your DADA2 seqtab output
taxonomy <- read_excel("~/Desktop/incubation_16S_final.xlsx", sheet = "tax_final (2)", na = c("","NA")) # fill this in with the tax_final file (DADA2 output but manipulated in Excel first)
metadata <- read_excel("~/Desktop/incubation_16S_final.xlsx", sheet = "metadata_final") # fill this in with the metadata file 

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

# Sub-setting the blanks that we want to remove (taxa)
dorm_subtract <- subset_samples(dorm1, sample_blank == "BLANK")

dorm_subtract

# Pseudo code 
# Tell R to subtract the blanks from the samples (in terms of ASV, remove ASV in the blanks from the samples)
# Vectorized approach (coding session with Hannah)
# step 1, convert to 0 and 1
mat <- otu_table(dorm_subtract)
mat <- as.matrix(data.frame(mat))
mat[mat > 0] <- 1
# step 2: come up with the row sum 
matrow <- rowSums(mat)
# remove anything from matrow that is 0

matrow_zero <- matrow[matrow == 0]
#sort(matrow_nozero)

goodOTUs <- names(matrow_zero)
dorm_taxaremoved <- prune_taxa(goodOTUs, dorm1)

hist(matrow_zero)

# Visualize the rarity curve
rarecurve(t(otu_table(dorm_taxaremoved)), step = 50, cex = 0.5, label = FALSE)

# Visualize via histogram
hist(sort(colSums(otu_table(dorm_taxaremoved)))) 

# rarefy without replacement
dorm1.rarefied_taxaremoved = rarefy_even_depth(dorm_taxaremoved, rngseed = 1, sample.size = 5000, replace = FALSE) ## removed sample sums but may put back

# Visualize the newly rarified curve
rarecurve(t(otu_table(dorm1.rarefied_taxaremoved)), step = 50, cex = 0.5, label = FALSE)

# Create the OTU table with the rarified data
total_otu_table <- otu_table(dorm1.rarefied_taxaremoved)

# Calculate relative abundance if appropriate
total = median(sample_sums(dorm1.rarefied_taxaremoved))
standf = function(x, t=total) round(t * (x / sum(x)))
dorm1 = transform_sample_counts(dorm_taxaremoved, standf)

# Preparing for ordination
dorm_transformed.taxaremoved <- t(sqrt(otu_table(dorm1.rarefied_taxaremoved)))## takes sq root of OTU table and does the t command, transform command (changes rows and columns)

# (Optional) Preparing for ordination 
##dorm_transformed <- dorm_transformed[rownames(dorm_transformed)!="insert any samples you want to pluck from your data here"

# Run a distance metric before ordination 
dm <- vegdist(dorm_transformed.taxaremoved, method = "bray")

# Run the ordinaton
dorm_taxaremoved.nmds <- metaMDS(dm, k = 2, trymax = 1500) ## k is the number of axes, sometimes 2 is not enough, trymax is the number of times it tries to find solutions
## stress values: how well the variation is represented (stress less than 0.05 is good, 0.3 is poor)

# Plot the MDS in vegan
plot(dorm_taxaremoved.nmds)
orditorp(dorm_taxaremoved.nmds,display = "sites",cex = 1.25,air = 0.01)

# Create a map data object
mapdata <- sample_data(dorm1.rarefied_taxaremoved)
mapdata$rowname <- rownames(mapdata)
dorm_taxaremoved.nmds.mapdata <- dorm_taxaremoved.nmds$points %>% data.frame() %>%
  rownames_to_column(var = "rowname") %>%   # Adds column called "rowname" to metadata and NMDS
  left_join(mapdata, by = "rowname") %>% # Matches row name from data to sample meta and joins them
  
 # select(sample_ID,
         #starts_with("NMDS1"), starts_with("NMDS2"), ## joined NMDS to mapping data by row name and re-ordered columns
         #MDS1, MDS2, everything())
# Making a "blank" column
#mutate(sample_blank = ifelse(test = grepl("BLANK", sample_ID), no = "Sample", yes = "BLANK")) %>%

# Get stress value
stress.dorm = paste("stress =", round(dorm_taxaremoved.nmds$stress, digits = 4)) ## stress = 0.1105

# Plotting the ordination in ggplot ## GETTING AN ERROR HERE ABOUT NO MAPDATA OBJECT 
ggplot(dorm_taxaremoved.nmds, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(),alpha = 1, size = 3) + # Alpha in "aes" shows more opaque (seethroughness)
  geom_text(aes(label = sample), size = 2)

# END

