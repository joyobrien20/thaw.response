# Joy O'Brien
# Microbial stats on Phyloseq object
# PRELIM
# April 29

# TO-DO: 
#   Make the top 5 phyla figure prettier 
#   NMDS / make this the hypoth 3 script (biplot and permanova)
#   Alpha diveristy 
# Load the necessary packages (this is from the phyloseq tutorial)
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("vegan")
library("viridis")
library("RColorBrewer")
library("colorRamps")
library("colorspace")



# FROM NATE'S 16S SCRIPT

phylum.sum = tapply(taxa_sums(dorm1rarefied), tax_table(dorm1rarefied)[, "Phylum"], sum, na.rm=TRUE)

top5phyla = names(sort(phylum.sum, TRUE))[1:5]
TP5 = prune_taxa((tax_table(dorm1rarefied)[, "Phylum"] %in% top5phyla), dorm1rarefied)

##taxa barplot
plot_bar(TP5, fill = "Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack")+
  scale_fill_manual(values = c("#ff9999", "#ffcc99", "#ffff99", "66b2ff", "99ffcc"))

# Make a bar graph of the data based on division ##PLOT TOP 5 PHYLA TO SEE THIS BETTER 
plot_bar(dorm1rarefied, fill = "Phylum")

# But make it prettier
plot_bar(dorm1rarefied, fill = "Phylum") + 
  geom_bar(aes(fill = Phylum, x = sample_ID), stat = "identity", position = "stack") +
  scale_fill_manual("Legend", values = c("Firmicutes" = "black"))

# Now let's make a basic heatmap
plot_heatmap(dorm1rarefied, method = "NMDS", distance = "bray")

### COME BACK TO THIS!

# Alpha diversity
plot_richness(dorm1rarefied, measures = c("Chao1", "Shannon"))
shannon <- estimate_richness(dorm1rarefied, measures = c("Chao1", "Shannon")) %>%
  geom_point(aes(size = 2, shape = 15))


## Alpha diversity take two (https://micca.readthedocs.io/en/latest/phyloseq.html)
# plot richness 
plot_richness(dorm1rarefied, x = "pre_post_thaw", color = "site", measures = c("Observed"))
# Make a boxplot of the number of OTUs and Shannon entropy 
plot_richness(dorm1rarefied, x = "site", measures = c("Observed", "Shannon")) + geom_boxplot()
estimate_richness(dorm1rarefied)


rich <- estimate_richness(dorm1rarefied)

#Stacked bar plot (https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)
dorm_phylum <- dorm1rarefied %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


# Plot with ggplot 
ggplot(dorm_phylum, aes(x = pre_post_thaw, y = Abundance, fill = Phylum)) + 
  facet_grid(site~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(
    breaks = c("7/8", "8/4", "9/2", "10/6"),
    labels = c("Jul", "Aug", "Sep", "Oct"), 
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") 
  #ggtitle("Phylum Composition of Lake Erie \n Bacterial Communities by Sampling Site") 


# Scale reads to even depth 
dorm_scale <- dorm1rarefied %>%
  scale_reads(round = "round") 

# Ordinate via PCOA
dorm_pcoa <- ordinate(
  physeq = dorm1rarefied, 
  method = "PCoA", 
  distance = "bray"
)

# Plot the PCOA
plot_ordination(
  physeq = dorm1rarefied,
  ordination = dorm_pcoa,
  color = "site",
  shape = "pre_post_thaw",
  title = "PCoA of permafrost bacterial communities pre- and post thaw"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


# Let's do an NMDS instead

# Run the NMDS and ordinate 
dorm_nmds <- ordinate(
  physeq = dorm1rarefied, 
  method = "NMDS", 
  distance = "bray"
)
print(dorm_nmds)
summary(dorm_nmds)
plot_ordination(
  physeq = dorm1rarefied,
  ordination = dorm_nmds,
  color = "site",
  shape = "pre_post_thaw",
  title = "NMDS of Permafrost Bacterial Communities"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

# BELOW IS AN EXAMPLE OF A PERMANOVA


# Calculate bray curtis distance matrix
dorm_bray <- phyloseq::distance(dorm1rarefied, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(dorm1rarefied))

# Adonis test
adonis(dorm_bray ~ site, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(dorm_bray, sampledf$site)
permutest(beta)

# Constrained ordination here 

## Alpha Diversity
min_lib <- min(sample_sums(dorm1rarefied))

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(dorm1rarefied)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(dorm1rarefied)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(dorm1rarefied)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(dorm1rarefied, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

alpha <- rbind(rich_stats, even_stats)
s <- data.frame(sample_data(dorm1rarefied))
alphadiv <- merge(alpha, s) 

ggplot(alphadiv, aes(x = site, y = mean, color = site, group = site, shape = site)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +
  facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) +
  scale_x_discrete(
    breaks = c("7/8", "8/4", "9/2", "10/6"),
    labels = c("Jul", "Aug", "Sep", "Oct"), 
    drop = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
