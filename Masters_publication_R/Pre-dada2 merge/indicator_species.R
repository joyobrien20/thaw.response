# Joy O'Brien
# Masters work in Ernakovich lab
# Attempting indicator species analysis

# Install if needed and then load the necessary library

library(indicspecies)
library(readxl)
library(phyloseq)


# Indicator species analysis 
# Load in csv ASV_transpose_v_2 (previously transposed in R and saved, and then downloaded for manipulating of sample IDs and removing DNA blanks)


# Defining the classification of sites (pre-thaw vs post-thaw)
groups = c(rep(1, 23), rep(2, 22))
groups

all_site_indic <- multipatt(ASV_transpose_v_2, groups, control = how(nperm=999))
summary(all_site_indic)

# Transposing the data frame
#dorm1 <-(otu_tablecsv)
#dorm1 <- as.matrix(otu_tablecsv)
#names(dorm1) <- as.matrix(dorm1[1, ]) # renames columns
#dorm1 <- dorm1[-1, ] #removes first row

#groups = c(rep(1, 17), rep(2, 14), rep (3,10))
#groups

#indval = multipatt(dorm1_num, groups, control=how(nperm=999))

#do some shuffling of the OTU table
dorm1 <- as.data.frame(t(otu_tablecsv)) #makes it a dataframe and puts x into y and y into x (flips it)

names(dorm1) <- as.matrix(dorm1[1, ]) # renames columns

dorm1 <- dorm1[-1, ] #removes first row
dorm1_num<-as.data.frame(lapply(dorm1, as.numeric)) #convert from character to number

dorm1_num$SampleID<-row.names(dorm1) #puts row names as sample ID column

#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
#metadata<-read.csv("Fixedmetadata2.csv") #read in metadata
#head(metadata) # check
# Read in metadata
meta <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final") # fill this in with the metadata file )
## Join based on SampleID
dorm1_Final<-left_join(dorm1, meta) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus = SpOTU_Final[,1:94] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat = SpOTU_Final$WaterType #the metadata column group you care about

SPind=multipatt(x=SPotus, cluster=SPwat, func = "r.g", control = how(nperm=9999))

summary(SPind)
