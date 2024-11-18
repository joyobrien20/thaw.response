# Joy O'Brien
# Master's project, Ernakovich lab
# Microbial stats on Phyloseq object
# May 7, 2023

# RE- sequenced updated script

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

# Read in the phyloseq object "dorm1rarefied_reseq"
# Read in the shannon diversity values that were previously obtained from this script but saved as an excel file and loaded back in 

shannon <- read_excel("~/GitHub/Masters_publication/shan_simp.xlsx", sheet = shannonsimp_reseq)

# Creating a table of Simpson and Shannon diversity
# THE FOLLOWING IS ONLY ON PRE THAW SAMPLES 

# Subset pre-thaw samples from the original phyloseq object (if needed)
presamples <- subset_samples(dorm1rarefied_reseq, pre_post_thaw == "pre")
pre_meta <- subset_samples(dorm1rarefied_reseq@sam_data, pre_post_thaw == "pre")

shansimp_pre <- estimate_richness(presamples, measures = c("Simpson", "Shannon"))
print(shansimp_pre)

# Subset samples post thaw for shannon and simpson diversity 
postsamples <- subset_samples(dorm1rarefied_reseq, pre_post_thaw == "post")
post_meta <- subset_samples(dorm1rarefied_reseq@sam_data, pre_post_thaw == "post")

# Run diversity analysis
shansimp_post <- estimate_richness(postsamples, measures = c("Simpson", "Shannon"))
print(shansimp_post)

# Shannon
shan_pre <- estimate_richness(presamples, measures = c("Shannon"))
print(shan_pre)

# Visualize pre-thaw shannon
plot_richness(dorm1rarefied_reseq, x = "site", color = "pre_post_thaw", measures = c("Simpson")) +
  geom_boxplot() +
  theme_classic() +
  xlab("Site") +
  theme(text = element_text(size =12)) +
  scale_color_discrete("Treatment")

# Let's see if it's normally distributed
hist(shansimp_pre$Shannon, main = "Shannon index", xlab = "")
shapiro.test(shansimp_pre$Shannon) # normally distributed woo

#***********************************************************************

# Working with distance matrices 
# Creating a bray curtis distance matrix from our phyloseq object
dist.bc <- distance(dorm1rarefied_reseq, method = "bray")
print(dist.bc)
# Unifrac
dist.uf <- distance(dorm1rarefied_reseq, method = "unifrac")
ord <- ordinate(dorm1rarefied_reseq, method = "MDS", distance = "bray")
# Visualize
plot_ordination(dorm1rarefied_reseq, ord, color = "site", shape = "pre_post_thaw") +
  stat_ellipse(aes(group = site))

#***************************************************************************************************
# Ordinating with environmental variables
 
# Getting rid of the column of numbers in the matrix
abs_abund <- abs_abund %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')
# convert to matrix?
abs_abund_mat <- as.matrix(abs_abund)
# Run the NMDS and ordinate 
dorm_nmds <- ordinate(
  physeq = dorm1rarefied_reseq, 
  method = "NMDS", 
  distance = "bray"
)
print(dorm_nmds)
summary(dorm_nmds)
env <- incubation_physical_chemical
en = envfit(dorm_nmds, env, permutations = 999, na.rm = TRUE)


# ATTEMPTING TO ADD IN SPECIES VECTORS
library(vegan)
library(ggplot2)
library(grid)
library(readxl)

abs_abund <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "abs_abund_all")

# Getting rid of the column of numbers in the matrix
abs_abund <- abs_abund %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')
# convert to matrix?
abs_abund_mat <- as.matrix(abs_abund)

# NMDS code for bi-plot? # COME BACK TO THIS AND MAKE IT WORK
nmds <- metaMDSdist(abs_abund_mat, distance = "bray")
nmds

env <- incubation_physical_chemical
en <- envfit(dorm_nmds,env, permutations = 999, na.rm = TRUE)
en
plot(nmds)

plot(en)

# NMDS (PUBLICATION FIGURE + labels)
plot_ordination(
  physeq = dorm1rarefied_reseq,
  ordination = dorm_nmds,
  color = "site",
  shape = "pre_post_thaw",
  #title = "NMDS of Permafrost Bacterial Communities"
) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(size = 2)+
  geom_text(aes(label = core), size = 4, color = "black") +
  theme_classic() +
  theme(axis.text = element_text(size = 18)) +
  scale_color_discrete("Site") +
  scale_shape_discrete("Treatment")
  scale_shape_manual("Treatment", labels = c("Pre-thaw", "Post-thaw"), values = c(2, 1)) +

# NMDS (PUBLICATION FIGURE + NO LABELS)
plot_ordination(
  physeq = dorm1rarefied_reseq,
  ordination = dorm_nmds,
  color = "site",
  shape = "pre_post_thaw",
  #title = "NMDS of Permafrost Bacterial Communities"
) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 18)) +
  scale_color_discrete("Site") +
  scale_shape_manual("Treatment", labels = c("Pre-thaw", "Post-thaw"), values = c(2, 1))

plot_ordination(
  physeq = dorm1rarefied_reseq,
  ordination = dorm_nmds,
  color = "site",
  shape = "pre_post_thaw",
) +
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 18)) +
  scale_color_discrete("Site", breaks = c("Site2", "Site1")) +
  scale_shape_discrete("Treatment", breaks = c("post-thaw", "pre-thaw"))

#****
dorm1rarefied_reseq$pre_post_thaw <- factor(dorm1rarefied_reseq$pre_post_thaw, levels = c("pre-thaw", "post-thaw"))
  print(dorm_nmds)
# colour = "grey90"
#**************************************************************************
# Calculate bray curtis distance matrix
dorm_bray <- phyloseq::distance(dorm1rarefied, method = "bray")

# Making a data frame from the phyloseq object 
dormOTU <- dorm1rarefied_OTU %>%
  data.frame()

# Saving matrix
# write.csv(dormOTU,"~/Desktop/dormOTU.csv")

# Need to transform so species are columns and samples are rows (for vegan)
dormOTUtransform <- t(sqrt(dormOTU))
dorm_bray <- vegdist(dormOTUtransform, method = "bray")

# Make a data frame from the sample_data
# sampledf <- data.frame(sample_data(dorm1rarefied))

# Adonis test (PERMANOVA)
# adonis(dorm_bray ~ pre_post_thaw, data = sampledf)

# Homogeneity of dispersion test for site
beta <- betadisper(dorm_bray, dorm1rarefied_sam_data$site) # testing the dispersion between groups (sites)
print(beta) # dispersion is how scattered the points are

adonis2(dorm_bray ~ site, data=datadorm)

# Homogeneity of dispersion test for core
beta_core <- betadisper(dorm_bray, dorm1rarefied_sam_data$core)
print(beta_core)

# Homogeneity of dispersion test for pre post thaw 
beta_pp <- betadisper(dorm_bray, dorm1rarefied_sam_data$pre_post_thaw)
print(beta_pp)
# Are these distances (dispersion/distance to the median) significant?

# Site 
permutest(beta)

# Core 
# permutest(beta_core)

# Pre post thaw 
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
# Making a data frame for PERMANOVA analysis
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

# Permanova assessing composition by core and pre post thaw 
corethaw.pnova <- vegan::adonis2(dorm_bray ~ core + pre_post_thaw, data = datadorm)
print(corethaw.pnova)

sitethaw.strata <- vegan::adonis2(dorm_bray ~ pre_post_thaw, permutations = perms, data = datadorm)
print(sitethaw.strata)

# SAVE
saveRDS(dorm1rarefied, "~/Desktop/dorm1.rds")
saveRDS(dorm_bray, "~/Desktop/dormbray.rds")

#Alpha Diversity (this code was not used for analysis)
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

#*******************************************************************

# Creating a table of Simpson and Shannon diversity and running stats

# Run diversity analysis
shansimp <- estimate_richness(dorm1rarefied_reseq, measures = c("Simpson", "Shannon"))
print(shansimp)

shan <- estimate_richness(dorm1rarefied, measures = c("Shannon"))
# Export as a CSV
# write.csv(shansimp,"~/Desktop/shansimp.csv")

# Reading all of the sites and samples and their shannon and simpson values into R
ssmeta <- read_excel("/Users/joyobrien/GitHub/Masters_publication/shan_simp.xlsx", sheet = "shansimp_reseq")


# Reading each site in individually to get significance codes pre and post thaw by site
sscrrel <- read_excel("/Users/joyobrien/GitHub/Masters_publication/shan_simp.xlsx", sheet = "shansimp_crrel")
ssfl <- read_excel("/Users/joyobrien/GitHub/Masters_publication/shan_simp.xlsx", sheet = "shansimp_fl")
ssut <- read_excel("/Users/joyobrien/GitHub/Masters_publication/shan_simp.xlsx", sheet = "shansimp_ut")


# Let's see if Shannon is normally distributed
hist(ssmeta$Shannon, main = "Shannon index", xlab = "")
shapiro.test(ssmeta$Shannon) # W = 0.89884, p-value = 0.001154, this is not normally distributed

# Testing for each site 
shapiro.test(sscrrel$Shannon) # p = 0.3176 normally distributed, use anova 
shapiro.test(ssfl$Shannon) # p = 0.1048 normally distributed, use anova 
shapiro.test(ssut$Shannon) # p = 0.0009 not normally distributed, use kruskal 

# Shannon pre-post-thaw by site 

# CRREL
crrelaov <- aov(Shannon ~ treatment, data = sscrrel)
summary(crrelaov)
crreltukey <- TukeyHSD(crrelaov, conf.level = 0.95)
print(crreltukey)

library(multcompView) #used to obtain the significance values
tukey_crrel <-multcompLetters4(crrelaov, crreltukey)
print(tukey_crrel)

# FL
flaov <- aov(Shannon ~ treatment, data = ssfl)
summary(flaov)
fltukey <- TukeyHSD(flaov, conf.level = 0.95)
print(fltukey)

tukey_fl <-multcompLetters4(flaov, fltukey)
print(tukey_fl)

# Utqiagvik
kruskal.test(ssut$Shannon ~ treatment, data = ssut)
kruskal(ssut$Shannon, ssut$treatment, group=TRUE, p.adj="bonferroni")$groups

# Obtain significance codes for Shannon and Simpson diversity

# Shannon by site
kruskal.test(Shannon ~ site, data = ssmeta)# chi-squared = 16.329, df = 2, p-value = 0.0002846
dunnTest(ssmeta$Shannon ~ site, data = ssmeta)

# Shannon by pre thaw 
kruskal.test(ssmeta$Shannon ~ thaw_type, data = ssmeta) # chi-squared = 14.702, df = 1, p-value = 0.0001259
kruskal(ssmeta$Shannon, ssmeta$thaw_type, group=TRUE, p.adj="bonferroni")$groups # to get significance codes

# Shannon by pre- post-thaw by site 
# CRREL
kruskal.test(ssmeta$Shannon ~ thaw_type, data = ssmeta) # chi-squared = 6.5455, df = 1, p-value = 0.01052
kruskal(ssmeta$Shannon, ssmeta$pre_post_thaw, site == "CRREL", group=TRUE, p.adj="bonferroni") #$groups

# Farmers Loop 
kruskal.test(ssmeta$Shannon ~ thaw_type, site == "FL", data = ssmeta) # chi-squared = 6, df = 1, p-value = 0.01431

# Utqiagvik
krusk.ut <- kruskal.test(ssmeta$Shannon ~ thaw_type, site == "Utqiagvik", data = ssmeta) # chi-squared = 8.0769, df = 1, p-value = 0.004483


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

# Shannon and Simpson diversity measurements visualization (Need to add in p-values)

shannonsimp_reseq <- read_excel("~/GitHub/Masters_publication/shan_simp.xlsx")
library(ggsignif)
# Plotting Shannon
ggplot(shannonsimp_reseq, aes(x = site, y = Shannon, color = thaw_type)) +
  geom_boxplot() + 
  #geom_signif(comparisons = list(c("Group1", "Group2")), annotations = "p=0.005") + trying to add significance to plot itself?
  theme_classic() +
  xlab("Site") +
  ylab("Alpha Diversity Measure") +
  scale_color_discrete("Treatment", labels=c('Pre-thaw', 'Post-thaw')) +
  theme(text = element_text(size = 12))

# Repeat above code with Simpson diversity
ggplot(shannonsimp_reseq, aes(x = site, y = Simpson, color = thaw_type)) +
  geom_boxplot() + 
  theme_classic() +
  xlab("Site") +
  ylab("Alpha Diversity Measure") +
  scale_color_discrete("Treatment", labels=c('Pre-thaw', 'Post-thaw')) +
  theme(text = element_text(size = 12))
#**********************************************************************************
# Extra code used for visualizing data

# Looking at the composition of the phyloseq object, we start with the top 5 phyla here: 
# Obtain the top 5 phyla
phylum.sum = tapply(taxa_sums(dorm1rarefied_reseq), tax_table(dorm1rarefied_reseq)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
TP5 = prune_taxa((tax_table(dorm1rarefied_reseq)[, "Phylum"] %in% top5phyla), dorm1rarefied_reseq)

# Visualize the top 5 phyla 
plot_bar(TP5, fill = "Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  theme_classic() 

# Now let's look at the top 10 phyla
# Obtain top 10 phyla 
top10phyla <- names(sort(phylum.sum, TRUE))[1:10]
TP10 <- prune_taxa((tax_table(dorm1rarefied_reseq)[, "Phylum"] %in% top10phyla), dorm1rarefied_reseq)

# Visualize top 10 phyla
plot_bar(TP10, fill = "Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  theme_classic()

# Obtain the top 5 Classes (for fun)
class.sum = tapply(taxa_sums(dorm1rarefied_reseq), tax_table(dorm1rarefied_reseq)[, "Class"], sum, na.rm=TRUE)
top5class = names(sort(class.sum, TRUE))[1:5]
CL5 = prune_taxa((tax_table(dorm1rarefied_reseq)[, "Class"] %in% top5class), dorm1rarefied_reseq)

# Visualize top 5 classes
plot_bar(CL5, fill = "Class") #+
#geom_bar(aes(fill="Genus"), stat="identity", position="stack")+
#scale_fill_manual(values = c("#ff9999", "#ffcc99", "#ffff99", "66b2ff", "99ffcc"))

# Now let's look at all of the phyla in the phyloseq object 
plot_bar(dorm1rarefied_reseq, fill = "Phylum")

#**************************************************************
#*# TOP 5 PHYLA IN PRE THAW SAMPLES 

# Obtain the top 5 phyla
phylum.sum_pre <- tapply(taxa_sums(presamples), tax_table(presamples)[, "Phylum"], sum, na.rm=TRUE)
top5phyla_pre = names(sort(phylum.sum, TRUE))[1:5]
TP5_pre = prune_taxa((tax_table(presamples)[, "Phylum"] %in% top5phyla), presamples)

# Let's visualize 
plot_bar(TP5_pre, fill = "Phylum") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  theme_classic()

# Top 10 phyla in pre-thaw samples 
phylum.sum_pre <- tapply(taxa_sums(presamples), tax_table(presamples)[, "Phylum"], sum, na.rm=TRUE)
top10phyla_pre = names(sort(phylum.sum_pre, TRUE))[1:10]
TP10_pre = prune_taxa((tax_table(presamples)[, "Phylum"] %in% top10phyla_pre), presamples)

# Visualize the top 10 phyla in pre-thaw samples by site
plot_bar(TP5_pre, x = "site", fill = "Phylum") +
  geom_bar(aes(color = Phylum , fill = Phylum), stat="identity", position = "stack") +
  theme_classic() +
  #scale_color_viridis(discrete = TRUE, option = "A") +
  display.brewer.all(n = 10, type = "all", select = NULL, colorblindFriendly = TRUE) +
  xlab("Site")

# Visualize the top 10 phyla in pre-thaw samples by core
plot_bar(TP10_pre, x = "site", fill = "Phylum") +
  geom_bar(aes(color = Phylum , fill = Phylum), stat="identity", position = "stack") +
  theme_classic() +
  #scale_color_viridis(discrete = TRUE, option = "A") +
  display.brewer.all(n = 10, type = "all", select = NULL, colorblindFriendly = TRUE) +
  xlab("Site") +
  ylab("Relative Abundance")+
  theme(text = element_text(size = 18))


# END
