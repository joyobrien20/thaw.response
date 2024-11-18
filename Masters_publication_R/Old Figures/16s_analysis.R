write('PATH="${RTOOLS40_HOME}(\\)usr(\\)bin;${PATH}"', file = "~/.Renviron", append = TRUE)

Sys.which("make")


install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics
install.packages("stringr")
install.packages("Biostrings")

## Load the packages
library("BiocManager")
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("Biostrings")
library("vegan")


str(RNA_meta)
str(OTU_RNA)



otu_mat <- OTU_RNA_DNA %>%
  tibble::column_to_rownames("OTU") 

tax_mat <- Taxonomy %>% 
  tibble::column_to_rownames("OTU")

samples_df <- Thesis_meta %>% 
  tibble::column_to_rownames("sample") 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

PO <- phyloseq(OTU, TAX, samples)
PO

#and we should get rid of all zero abundance samples as they will mess up everything

PO1 <- prune_samples(sample_sums(PO)>0,PO)
PO1


#There are some taxa that are unassigned past the Kingdom level.  These are removed from the analysis.  

PO1 <- subset_taxa(PO1, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
PO1

###**Analysis note: Decision was made on 21st November that all chloroplasts distinctive from cyanobacteria will be removed from all analysis and datasets. Chloroplasts are distinguished from cyanobacteria at the class level


PO1 <- subset_taxa(PO1, Class!="NA"| is.na(Class))


#Remove OTUs that do not show appear more than 5 times in more than half the samples

wh0 = genefilter_sample(ps1, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps1))
A = prune_taxa(wh0, ps1)

##relative abundance if appropriate
total = median(sample_sums(PO1))
standf = function(x, t=total) round(t * (x / sum(x)))
PO1 = transform_sample_counts(PO, standf)

##n to rarify, vegan package
rarecurve(t(otu_table(PO1)), step=50, cex=0.5)

# rarefy without replacement
PO1.rarefied = rarefy_even_depth(PO1, rngseed=1, sample.size=1*min(sample_sums(PO1)), replace=F)


#Keep only the most abundant five phyla.

phylum.sum = tapply(taxa_sums(PO1.rarefied), tax_table(PO1.rarefied)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
TP5 = prune_taxa((tax_table(PO1.rarefied)[, "Phylum"] %in% top5phyla), PO1.rarefied)

##taxa barplot
plot_bar(TP5, fill = "Phylum")+
  geom_bar(aes(fill=Phylum), stat="identity", position="stack")+
  scale_fill_manual(values = c("#ff9999", "#ffcc99", "#ffff99", "66b2ff", "99ffcc"))

  

##alpha diversity
plot_richness(PO1.rarefied, measures=c("Chao1", "Shannon"))
shan <- estimate_richness(FOX, measures=c("Chao1", "Shannon"))

+ geom_point(size=2, shape=15)

PO1.ord <- ordinate(PO1.rarefied, "PCoA", "bray")
plot_ordination(PO1.rarefied, PO1.ord, type="sample", color="Site", shape = "DNA_RNA", label = "sample_ID" , 
                title="OTUs")
PO1.ord$values
FOX<-subset_samples(PO1,Core =="FOX1"|Core=="FOX2"|Core=="FOX3"|Core=="FOX4"|Core=="Blank")

FOX.ord <- ordinate(FOX, "PCoA", "bray")
plot_ordination(FOX, FOX.ord, type="sample", color="Core", shape = "DNA_RNA", label = "sample_ID",
                title="OTUs") 

metadata <- as(sample_data(FOX), "data.frame")

adonis(distance(FOX, method="bray") ~ Core,
       data = metadata)


#fit regression model
model <- lm(Euler~ShannonD, data=R_Format.3)
summary(model)

a <-ggplot(R_Format.3, aes(x=Euler, y=ShannonD))+
  geom_point(size=2, shape=15, color = "red")+
  geom_smooth(method=lm, se=FALSE)+
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.text.y = element_text(size=15)) +
  theme(axis.title = element_text(size=15))
a + ylab("shannon")
a + xlab("%")+ ylab("shannon")

test <- otu_table(PO1)
test

sort(colSums(test))

## transform phyloseq otu_table
sb_transformed <- t(sqrt(test+.00001))
dm_bc <- vegdist(sb_transformed, method = "bray")
dm_bc

##mantel test
bc_uf_comp.mantel <- mantel(dm_bc, dm_uf, method = "pearson")

##log2fold change ##looking at differences between RNA and DNA communities
##multipatt indicspecies

##need to rarify, vegan package
rarecurve(t(otu_table(PO1)), step=50, cex=0.5)

# rarefy without replacement
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)



### Indicator species
install.packages("indicspecies")
library(indicspecies)

abund = OTU_RNA_DNA_indic_species[,5:ncol(OTU_RNA_DNA_indic_species)]
DNA_RNA = OTU_RNA_DNA_indic_species$DNA_RNA
Site = OTU_RNA_DNA_indic_species$Site
Core = OTU_RNA_DNA_indic_species$Core


inv = multipatt(abund, DNA_RNA, func = "r.g", control = how(nperm=9999))
summary(inv)

##Filter out indicator species from output
##taxonomy_table %>% filter(ASV %in% column_of_ASVs)

taxonomy_table_indicDNA <- Taxonomy %>% filter(OTU %in% indic_DNA$OTU)
taxonomy_table_indicDNA

library(tidyverse)
install.packages(tidyverse)
indic_DNA <- read_csv("~/Downloads/multipattDNA.csv")

Taxonomy <- read_csv("~/Downloads/Taxonomy.csv")

taxonomy_table_indicDNA <- Taxonomy %>% 
  filter(OTU %in% multipattDNA$OTU)

dim(taxonomy_table_indicDNA)
dim(Taxonomy)

taxonomy_table_indicDNA



##Top10
b1_genus<-tax_glom(PO1.rarefied, taxrank="Genus")
Genus10 = names(sort(taxa_sums(b1_genus), TRUE)[1:10])

##taxa barplot
plot_bar(Genus10, fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
