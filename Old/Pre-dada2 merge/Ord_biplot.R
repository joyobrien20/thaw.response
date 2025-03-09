# Ordinating with vectors
https://david-barnett.github.io/microViz/articles/web-only/ordination.html 


library(phyloseq)
library(ggplot2)
library(microViz)
library(vegan)
# install.packages("microViz")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
# install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)
# Load phyloseq object 
dorm <- dorm1rarefied_reseq
# Fix any uninformative taxonomic entries
dorm <- tax_fix(dorm) # try tax_fix_interactive if you have problems with your own data
ibd <- phyloseq_validate(dorm, remove_undetected = TRUE)

# Consider transforming microbial counts (using center log rank transformation which is recommended for seq. data)
ibd %>%
  tax_transform(trans = "clr", rank = "Class")

# PCoA requires a matrix of pairwise distances between samples 
# Normally you should NOT transform your data when using a distance-based method, but 
# it is useful to record an “identity” transformation anyway, to make it clear you have not transformed your data.
ibd %>%
  tax_transform(trans = "identity", rank = "Class") %>%
  dist_calc("bray")

ibd %>%
  tax_transform("clr", rank = "Class") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "pre_post_thaw", shape = "site", size = 2) +
  scale_colour_brewer(palette = "Dark2")


ibd %>%
  tax_transform("clr", rank = "Class") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCA") %>%
  ord_plot(color = "pre_post_thaw", shape = "site", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "Dark2")

# Code used in Berhardt et al., 2022
# https://online.ucpress.edu/elementa/article/10/1/00023/194712/Soil-microbial-communities-vary-in-composition-and
###-------------------------Community Analysis---------------------------###

####SHANNONS DIVERSITY#### - FT, NT and ALL
```{r}
#checking out shannon's diversity for class and species between no till, full-till seperately and the whole data set

shan_cl_NT <- diversity(jst_cl_NT)
species_cl_NT$cl_shan <- shan_cl_NT
NTshan.anov <- aov(cl_shan ~ AggSize, data = species_cl_NT)
summary(NTshan.anov)


shan_cl_FT <- diversity(jst_cl_FT)
species_cl_FT$cl_shan <- shan_cl_FT
FTshan.anov <- aov(cl_shan ~ Pesticide, data = species_cl_FT)
summary(FTshan.anov)

shan_cl <- diversity(jst_cl)
species_cl$cl_shan <- shan_cl
SPshan.anov <- aov(cl_shan ~ AggSize * Tillage, data = species_cl)
summary(SPshan.anov)
#Shannon diversity is significantly different between aggregate size classes when full-till and no-till plots are combined together into one analysis
```


####NMDS/PERMANOVA#### - ALL - class

#log transforming my bacterial community data and making a distance matrix from the results
#log+1 used because log(0) is undefined 
```{r}
lg_jst_cl <- log1p(jst_cl)
lg_jst_cl <- data.matrix(lg_jst_cl)

cl_bray_dist <- vegdist(jst_cl, method = "bray")

```
# Joy's attempt: 
dorm_mat <- as.matrix(dorm1rarefied_OTU_reseq)
lg_dorm_mat <- log1p(dorm_mat)
lg_dorm_mat <- data.matrix(lg_dorm_mat)

bray_dist <- vegdist(lg_dorm_mat, method = "bray")
nmds <- metaMDS(bray_dist)

#NMDS: with bray curtis distance measure. Performed for multiple axes. axis where stress decreased by >0.05 was retained. Best solution is on three axes
Stress evalution:
  <.05 - excellent
.05-.10 - good
.10-.20 - usable, potential to mislead at upper end

#Running NMDS at with different amount of axes to see which one has the best stress (bray curtis distance used)
```{r}
NMDS_lg_cl_1 <- metaMDS(lg_dorm_mat, distance = "bray", k= 1, autotransform = FALSE, trymax = 1000)
NMDS_lg_cl_2 <- metaMDS(lg_dorm_mat, distance = "bray", k=2, autotransform = FALSE, trymax= 1000)
NMDS_lg_cl_3 <- metaMDS(lg_dorm_mat, distance = "bray", k=3, autotransform = FALSE, trymax = 1000)
NMDS_lg_cl_4 <- metaMDS(lg_dorm_mat, distance = "bray", k= 4, autotransform = FALSE, trymax = 1000)
NMDS_lg_cl_5 <- metaMDS(lg_dorm_mat, distance = "bray", k=5, autotransform = FALSE, trymax = 1000)
NMDS_lg_cl_6 <- metaMDS(lg_dorm_mat, distance = "bray", k=6, autotransform = FALSE, trymax = 1000)


#determined 3 axis solution was best
spp_stressvector<-as.vector(c(NMDS_lg_cl_1$stress, NMDS_lg_cl_2$stress, NMDS_lg_cl_3$stress, NMDS_lg_cl_4$stress, NMDS_lg_cl_5$stress, NMDS_lg_cl_6$stress))
plot(spp_stressvector)

#looking at the stress for solutions with axes 2-5
NMDS_lg_cl_2$stress
NMDS_lg_cl_3$stress
NMDS_lg_cl_4$stress
NMDS_lg_cl_5$stress

```
***Adding environmental variables***
  
  #Making environmental matrix should just include the metadata (block, plot etc...) and any other environmental variables of interest (pH, EC, total carbon, etc...)
  ```{r}
head(species_cl)
env_var <- as.data.frame(species_cl[,1:11])
env_var
```


#extracting the site and species scores from the NMDS for plotting
#Choose which NMDS axis you want to plot and substitute for "NMDS_lg_cl_3"
```{r}
site.scrs <- as.data.frame(scores(NMDS_lg_cl_3, display = "sites"))

#Two ways to accomplish the same thing

site.scrs$AggSize = env_var$AggSize
site.scrs$Tillage = env_var$Tillage
site.scrs$Block = env_var$Block
site.scrs

#Species Scores
spp.scrs <- as.data.frame(scores(NMDS_lg_cl_3, display = "species"))
spp.scrs <- cbind(spp.scrs, Species = row.names(spp.scrs))

spp.scrs
```

```{r}
### First, let's relativize the environmental data. We only actually have one numeric environmental variables, so it is not important in this situation, but when you have several environmental variables, it frequently is. Imagine environmental variables whose units are not comparable (ie. pH, depth of soil horizon, canopy cover, etc.). It is particularly important to standardize such values. The vegan function decostand() is really useful for this purpose! 

env_var

env_var$Block <- as.factor(env_var$Block)
env_var$AggSize <- as.factor(env_var$AggSize)
env_var$Tillage <- as.factor(env_var$Tillage)
env_var$Pesticide <- as.factor(env_var$Pesticide)


env_var_rel_lg <- decostand(select_if(env_var, is.numeric), method = "log")

env_var_rel_fit <- envfit(NMDS_lg_cl_3, env_var_rel_lg, choices = 1:2, permutations = 999)
env.scores <- as.data.frame(scores(env_var_rel_fit, display = "vectors"))
env.scores
env.scores2 <- cbind(env.scores, env.variables = rownames(env.scores), pval = env_var_rel_fit$vectors$pvals)
env.scores2 <- subset(env.scores2, pval<=0.05)
env.scores2
# NMDS_1v2_spp_scrs_with_env
```
***Plotting***
  
  Add sites
```{r}
plot(NMDS_lg_cl_3, "sites")
orditorp(plot(NMDS_lg_cl_3), display="sites", col= "red", air= 0.01)   # Gives points labels

#color by till_treat
NMDS_1v2 <- ggplot(site.scrs, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(NMDS1, NMDS2, shape = Tillage, colour = factor(AggSize)), size = 2.6) + labs(colour = "Aggregate Size Class", shape = "Tillage Treatment") + theme_classic() + theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 8, face = "bold")) + scale_color_manual(values= c("#ffcc00", "#669900", "#663300", "#e60000"))
NMDS_1v2

species_cl


```


#Add environmental variables
```{r}
NMDS_1v2_with_env <- NMDS_1v2 + geom_segment(data = env.scores2, aes(x = 0, xend = NMDS1*0.85, y = 0, yend = NMDS2*0.85), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd = 0.3) + ggrepel::geom_text_repel(data = env.scores2, aes(x = NMDS1*0.9, y = NMDS2*0.9, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)
NMDS_1v2_with_env                                                                  
ggsave("16S_community_NMDS.png", plot = NMDS_1v2_with_env, dpi = 400, height = 5, width = 7)



```



#Saving figures for publication!
```{r}
NMDS_1v2_with_env
setwd("_______________________")
ggsave("community_class_NMDS.png", plot = NMDS_1v2_with_env, dpi = 400, height = 5, width = 6)

```


By changing NMDS1 and NMDS2, you can make figures for NMDS1 vs. NMDS3 and NMDS2 vs NMDS3

NMDS axis: 
  The interpretation of NMDS ordinations is fundementally different than PCA. The first axis in PCA, by definition, always explains the most variation in the data. The second axis explains the second most variation orthoganal to the first axis. NMDS axis should be treated as meaningless. The number of axis is defined by the user, and the algorithm then attempts to reposition the observations in a way that maintains their ranked distances.  


********************************PERMANOVA*********************************************************
  
  ```{r}
#Adonis acts on matrices
lg_cl_matrix <- as.matrix(jst_cl)


