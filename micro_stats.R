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
library("ggpubr")
library("readxl")
library("FSA")
library("agricolae") # used to get significance codes 
# Looking at the phyloseq object

# Obtain the top 5 phyla
phylum.sum = tapply(taxa_sums(dorm1rarefied), tax_table(dorm1rarefied)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
TP5 = prune_taxa((tax_table(dorm1rarefied)[, "Phylum"] %in% top5phyla), dorm1rarefied)

# Obtain the top 5 class
class.sum = tapply(taxa_sums(dorm1rarefied), tax_table(dorm1rarefied)[, "Class"], sum, na.rm=TRUE)
top5class = names(sort(class.sum, TRUE))[1:5]
CL5 = prune_taxa((tax_table(dorm1rarefied)[, "Class"] %in% top5class), dorm1rarefied)

# Visualize top class
plot_bar(CL5, fill = "Class") #+
  #geom_bar(aes(fill="Genus"), stat="identity", position="stack")+
  #scale_fill_manual(values = c("#ff9999", "#ffcc99", "#ffff99", "66b2ff", "99ffcc"))
                    
# Visualize with a taxa barplot (Nate's script)
plot_bar(TP5, fill = "Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack")+
  scale_fill_manual(values = c("#ff9999", "#ffcc99", "#ffff99", "66b2ff", "99ffcc"))

# Looking at ALL of the phyla within the data
plot_bar(dorm1rarefied, fill = "Phylum")

# Cleaning up the phyla barplot
plot_bar(dorm1rarefied, fill = "Phylum") + 
  geom_bar(aes(fill = Phylum, x = sample_ID), stat = "identity", position = "stack") +
  scale_fill_manual("Legend", values = c("Firmicutes" = "black"))

#***********************************************************************
# Creating a table of Simpson and Shannon diversity
# THE FOLLOWING IS ONLY ON PRE THAW SAMPLES 
# Subset pre-thaw samples from the original phyloseq object (if needed)
presamples <- subset_samples(dorm1rarefied, pre_post_thaw == "pre")
pre_meta <- subset_samples(dorm1rarefied_sam_data, pre_post_thaw == "pre")

# Run diversity analysis
shansimp_pre <- estimate_richness(presamples, measures = c("Simpson", "Shannon"))
print(shansimp_pre)

# Visualize
boxplot(shansimp_pre)

# Let's see if it's normally distributed
hist(shansimp_pre$Shannon, main = "Shannon index", xlab = "")
shapiro.test(shansimp_pre$Shannon) # normally distributed woo

#***********************************************************************
# TOP 5 PHYLA IN PRE THAW SAMPLES 

# NOTE: if you end up using this figure in thesis/manuscript you need to make it cleaner
# Obtain the top 5 phyla
phylum.sum_pre <- tapply(taxa_sums(presamples), tax_table(presamples)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_pre = names(sort(phylum.sum, TRUE))[1:5]
TP5_pre = prune_taxa((tax_table(presamples)[, "Phylum"] %in% top5phyla), presamples)

# Visualize the top 5 phyla in pre-thaw samples
plot_bar(TP5_pre, x = "sample_ID", fill = "Phylum") +
  geom_bar(aes(color = Phylum , fill = Phylum), stat="identity", position = "stack")


#************************************************************************************
# TOP 5 PHYLA IN POST THAW SAMPLES 

# Subset post-thaw samples from the original phyloseq object
# NOTE: if you end up using this figure in thesis/manuscript you need to make it cleaner
postsamples <- subset_samples(dorm1rarefied, pre_post_thaw == "post")
shansimp_post <- estimate_richness(postsamples, measures = c("Simpson", "Shannon"))
print(shansimp_post)

# Obtain the top 5 phyla
phylum.sum_post <- tapply(taxa_sums(postsamples), tax_table(postsamples)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_post = names(sort(phylum.sum, TRUE))[1:5]
TP5_post = prune_taxa((tax_table(postsamples)[, "Phylum"] %in% top5phyla), postsamples)

# Visualize the top 5 phyla in post-thaw samples
plot_bar(TP5_post, x = "sample_ID", fill = "Phylum") +
  geom_bar(aes(color = Phylum , fill = Phylum), stat="identity", position = "stack") +
  theme_classic()
#************************************************************************************

# Working with distance matrices 
dist.bc <- distance(dorm1rarefied, method = "bray")
print(dist.bc)

dist.uf <- distance(dorm1rarefied, method = "unifrac")
ord <- ordinate(dorm1rarefied, method = "MDS", distance = "bray")

plot_ordination(dorm1rarefied, ord, color = "site", shape = "pre_post_thaw") +
  stat_ellipse(aes(group = site))

# Insert permanova code here? 




# Extra code for making a bar plot
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

# Ordination 

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
  #title = "NMDS of Permafrost Bacterial Communities"
) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)+
  theme_classic() +
  theme(text = element_text(size = 12)) +
  scale_color_discrete("Site") +
  scale_shape_discrete("Treatment")

print(dorm_nmds)
  
#**************************************************************************
# Calculate bray curtis distance matrix
dorm_bray <- phyloseq::distance(dorm1rarefied, method = "bray")

# Making a data frame from the phyloseq object 
dormOTU <- dorm1rarefied_OTU %>%
  data.frame()

# write.csv(dormOTU,"~/Desktop/dormOTU.csv")
# Need to transform so species are columns and samples are rows (for vegan)
dormOTUtransform <- t(sqrt(dormOTU))
dorm_bray <- vegdist(dormOTUtransform, method = "bray")

# make a data frame from the sample_data
# sampledf <- data.frame(sample_data(dorm1rarefied))

# Adonis test
# adonis(dorm_bray ~ pre_post_thaw, data = sampledf)

# Homogeneity of dispersion test for site
beta <- betadisper(dorm_bray, dorm1rarefied_sam_data$site) # testing the dispersion between groups (sites)
print(beta) # dispersion is how scattered the points are

adonis2(dorm_bray ~ site, data=datadorm)

# Homogeneity of dispersion test for core
beta_core <- betadisper(dorm_bray, dorm1rarefied_sam_data$core)
print(beta_core)

# Homogenity of dispersion test for pre post thaw 
beta_pp <- betadisper(dorm_bray, dorm1rarefied_sam_data$pre_post_thaw)
print(beta_pp)
# Are these distances (dispersion/distance to the median) significant?

# site 
permutest(beta)

# core 
# permutest(beta_core)

# pre post thaw 
permutest(beta_pp)

# Tukey for site 
tukey_beta <- TukeyHSD(beta)
plot(tukey_beta) 
print(tukey_beta)

# Tukey for core 
# tukey_beta_core <- TukeyHSD(beta_core)
# print(tukey_beta_core)

# Tukey for pre post thaw  # not necessary since it wasnt significant
tukey_beta_pp <- TukeyHSD(beta_pp)
print(tukey_beta_pp)

#****************************************************************************
# Making a data frame for permanova analysis
datadorm <- data.frame(sample_data(dorm1rarefied)) 

# site 
perms <- with(datadorm, how(nperm = 1000, blocks = site))
# trying this with core 
# perms_core <- with(datadorm, how(nperm = 1000, blocks = core))

# pre post thaw? 
perms_pp <- with(datadorm, how(nperm = 1000, blocks = pre_post_thaw))

# Permanova assessing composition by site and pre-post thaw 
sitethaw.pnova <- vegan::adonis2(dorm_bray ~ site + pre_post_thaw, data = datadorm)
print(sitethaw.pnova)

sitethaw.strata <- vegan::adonis2(dorm_bray ~ pre_post_thaw, permutations = perms, data = datadorm)
print(sitethaw.strata)

saveRDS(dorm1rarefied, "~/Desktop/dorm1.rds")
saveRDS(dorm_bray, "~/Desktop/dormbray.rds")
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

# Untouched code below from online 
# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Untouched code below from online 
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

#*******************************************************************
# Creating a table of Simpson and Shannon diversity and running stats

# Run diversity analysis
shansimp <- estimate_richness(dorm1rarefied, measures = c("Simpson", "Shannon"))
print(shansimp)

# Export as a CSV
# write.csv(shansimp,"~/Desktop/shansimp.csv")

# Reading all of the sites and samples and their shannon and simpson values into R
ssmeta <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final_shansimp")

# Reading each site in individually to get significance codes pre and post thaw by site 
sscrrel <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "shansimp_crrel")
ssfl <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "shansimp_fl")
ssut <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "shansimp_ut")

# Visualize
boxplot(shansimp)

# Let's see if Shannon is normally distributed
hist(ssmeta$Shannon, main = "Shannon index", xlab = "")
shapiro.test(ssmeta$Shannon) # W = 0.87838, p-value = 0.002171 this is not normally distributed

# Testing for each site 
shapiro.test(sscrrel$Shannon) # p = 0.2624 normally distributed, use anova 
shapiro.test(ssfl$Shannon) # p = 0.2112 normally distributed, use anova 
shapiro.test(ssut$Shannon) # p = 0.007 not normally distributed, use kruskal 

# Shannon pre-post-thaw by site 
# CRREL
crrelaov <- aov(Shannon ~ pre_post_thaw, data = sscrrel)
summary(crrelaov)
crreltukey <- TukeyHSD(crrelaov, conf.level = 0.95)
print(crreltukey)

library(multcompView)
tukey_crrel <-multcompLetters4(crrelaov, crreltukey)
print(tukey_crrel)
# FL

flaov <- aov(Shannon ~ pre_post_thaw, data = ssfl)
summary(flaov)
fltukey <- TukeyHSD(flaov, conf.level = 0.95)
print(fltukey)

tukey_fl <-multcompLetters4(flaov, fltukey)
print(tukey_fl)
# Utqiagvik 
kruskal.test(ssut$Shannon ~ pre_post_thaw, data = ssut)
kruskal(ssut$Shannon, ssut$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups

# Shannon by site
kruskal.test(ssmeta$Shannon ~ site, data = ssmeta)# chi-squared = 13.765, df = 2, p-value = 0.001026
kruskal(ssmeta$Shannon, ssmeta$site, group=TRUE, p.adj="bonferroni")$groups #using this to get the letter significance codes 
dunnTest(ssmeta$Shannon ~ site, data = ssmeta)

# Shannon by pre thaw 
kruskal.test(ssmeta$Shannon ~ pre_post_thaw, data = ssmeta) # chi-squared = 14.702, df = 1, p-value = 0.0001259
kruskal(ssmeta$Shannon, ssmeta$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups # to get significance codes

# Shannon by pre- post-thaw by site 
# CRREL
kruskal.test(ssmeta$Shannon ~ pre_post_thaw, site == "CRREL", data = ssmeta) # chi-squared = 6.5455, df = 1, p-value = 0.01052
kruskal(ssmeta$Shannon, ssmeta$pre_post_thaw, site == "CRREL", group=TRUE, p.adj="bonferroni") #$groups

# Farmers Loop 
kruskal.test(ssmeta$Shannon ~ pre_post_thaw, site == "FL", data = ssmeta) # chi-squared = 6, df = 1, p-value = 0.01431

# Utqiagvik
kruskal.test(ssmeta$Shannon ~ pre_post_thaw, site == "Utqiagvik", data = ssmeta) # chi-squared = 8.0769, df = 1, p-value = 0.004483

# Subsetting by site manually 
crrel_ss <- subset(ssmeta,subset =  site == "CRREL")
kruskal.test(crrel_ss$Shannon ~ pre_post_thaw, data = crrel_ss)

#dunnTest(ssmeta$Shannon ~ pre_post_thaw, data = ssmeta) 
#**********************************************************************************
# Let's see if Simpson is normally distributed
hist(ssmeta$Simpson, main = "Simpson index", xlab = "")
shapiro.test(shansimp$Simpson) # p-value = 1.541e-07, this is not normally distributed

# Checking to see if simpson is normally distributed when we work with each site individually 

shapiro.test(sscrrel$Simpson) # p = 0.01377 # use kruskal to evaluate 
shapiro.test(ssfl$Simpson) # p = 0.02356 # use kruskal
shapiro.test(ssut$Simpson) # p = 0.000399 not normally distributed, use kruskal 

# Simpson pre post thaw by site (obtaining the significance codes)

# CRREL
kruskal.test(sscrrel$Simpson ~ pre_post_thaw, data = sscrrel)
kruskal(sscrrel$Simpson, sscrrel$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups
# FL
kruskal.test(ssfl$Simpson ~ pre_post_thaw, data = ssfl)
kruskal(ssfl$Simpson, ssfl$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups
# Utqiagvik 
kruskal.test(ssut$Simpson ~ pre_post_thaw, data = ssut)
kruskal(ssut$Simpson, ssut$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups
# Simpson by site 
kruskal.test(ssmeta$Simpson ~ site, data = ssmeta) # chi-squared = 9.5593, df = 2, p-value = 0.008399
dunnTest(ssmeta$Simpson ~ site, data = ssmeta)

# Using the base Kruskal test we find that it's the same as above with the FSA package Kruskal test 
#stats::kruskal.test(ssmeta$Simpson ~ site, data = ssmeta)
# Checking to see if Simpson is normally distributed
shapiro.test(ssmeta$Simpson) # this is not normally distributed, p = 1.541e-07

# Simpson by site 
kruskal.test(ssmeta$Simpson ~ site, data = ssmeta)
kruskal(ssmeta$Simpson, ssmeta$site, group=TRUE, p.adj="bonferroni")$groups
dunnTest(ssmeta$Simpson ~ site, data = ssmeta)
# Simpson by pre- post-thaw 
kruskal.test(ssmeta$Simpson ~ pre_post_thaw, data = ssmeta) #chi-squared = 18.564, df = 1, p-value = 1.643e-05
kruskal(ssmeta$Simpson, ssmeta$pre_post_thaw, group=TRUE, p.adj="bonferroni")$groups

# Simpson by pre-post thaw by site 
# CRREL
kruskal.test(ssmeta$Simpson ~ pre_post_thaw, site == "CRREL", data = ssmeta) # chi-squared = 6.5455, df = 1, p-value = 0.01052

# FL
kruskal.test(ssmeta$Simpson ~ pre_post_thaw, site == "FL", data = ssmeta) # chi-squared = 6, df = 1, p-value = 0.01431

# Utqiagvik
kruskal.test(ssmeta$Simpson ~ pre_post_thaw, site == "Utqiagvik", data = ssmeta) # chi-squared = 8.0769, df = 1, p-value = 0.004483

#**********************************************************************************

#Alpha diversity take two (https://micca.readthedocs.io/en/latest/phyloseq.html)
# Alpha diversity: plotting Shannon and Chao1
plot_richness(dorm1rarefied, measures = c("Chao1", "Shannon"))

# Plot Shannon diversity/making the plot prettier
plot_richness(dorm1rarefied, x = "site", color = "pre_post_thaw", measures = c("Shannon")) +
  geom_boxplot(aes(x = reorder(pre_post_thaw))) + # need to reorder so that pre comes before post 
  theme_classic() +
  xlab("Site") +
  theme(text = element_text(size =12)) +
  scale_color_discrete("Treatment")


# Plot Simpson Diversity
plot_richness(dorm1rarefied, x = "site", color = "pre_post_thaw", measures = c("Simpson")) +
  geom_boxplot() +
  theme_classic() +
  xlab("Site") +
  theme(text = element_text(size =12)) +
  scale_color_discrete("Treatment")


plot_richness(dorm1rarefied, x = "site", color = "pre_post_thaw", measures = c("Simpson")) + 
  geom_boxplot()

# Make a boxplot of the number of OTUs and Shannon entropy 
plot_richness(dorm1rarefied, x = "site", measures = c("Observed", "Shannon")) + geom_boxplot()


