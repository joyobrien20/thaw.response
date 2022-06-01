# Joy O'Brien
# May 30, 2022
# Master's project
# Thaw response code (response ratio analysis code from Hannah Holland-Moritz)

# Load the necessary libraries

library(mctoolsr)
library(ggplot2)
library(tibble)
library(dplyr)
library(multcompView)
library(phyloseq)
library(readxl)
library(biomformat)
library(readr)

# Installing and using mctools, it wasn't as straight forward as install.packages
# install.packages("devtools")
# devtools::install_github("leffj/mctoolsr")

# Load rarified data and sample list 
#tax_table <- read_delim("~/Desktop/JMO_masters_2022/taxtable.txt")
#metadata <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final")
#OTU <-read_delim("~/Desktop/JMO_masters_2022/otu_table.txt")
#thaw.stats <- load_taxa_table("~/Desktop/JMO_masters_2022/taxtable.txt", "~/Desktop/JMO_masters_2022/otu_table.txt") # I cannot get this to work!
#mapping_file <- thaw.stats$map_loaded 

# Response ratio analysis ## need to edit with your own specs and file names! 
# Okay apparently I don't need to fill in my data in the above chunk-I can just load the table in here from phyloseq (-':

# read in file (load ASV table) 
thaw_response <- dormOTU
thaw_response<- thaw_response + 0.01
# Replace all 0 with NA
#thaw_response[thaw_response == 0] = 0.01
# creating a new variable (TOF) new table w/ 0 = NA ###
ToF <- thaw_response;
ToF[ToF > 0] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using relative abundance table, creates average
# treatment_names <- c("sample1", "sample2"...)
ave1 = apply(thaw_response[,c(1:12,24,30:31)],1,mean,na.rm = TRUE);# Control 
ave2 = apply(thaw_response[,c(13:23,25:29)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
sd1 = apply(thaw_response[,c(1:12,24,30:31)],1,sd,na.rm = TRUE);# Control #
sd2 = apply(thaw_response[,c(13:23,25:29)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
num1 = apply(ToF[,c(1:12,24,30:31)],1,sum,na.rm = TRUE); #Control#
num2 = apply(ToF[,c(13:23,25:29)],1,sum,na.rm = TRUE) # Treatment #

# Response ration calculations########################
# sample varriance #
v1 = sd1^2/(ave1^2*num1); # Control #
v2 = sd2^2/(ave2^2*num2) # treatment #

# square root of varriance # 
sqv = sqrt(v1+v2)

# RR = Response ratio #
RR = log(ave2/ave1)

# Confidence interval # 
ICl = RR-1.96*sqv
ICu = RR+1.96*sqv

# error #
err = ICu-RR
# Significance #
Sig = ICl*ICu
# creating data frame w/ 3 coloums (RR,err,Sig) #
res = as.data.frame(cbind(RR,err,Sig))
# create a new column from row names #
res$ID <- rownames(res)
# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# Come back to this and decide what to do when Hannah gets back to you
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"


write.table(res,"thaw_response.txt",sep = "\t")

tax_response <- left_join(x = res, y = tax_tablecsv, by = c("ID"="#ASV_ID"))
# save this 

# Filter tax response so that we don't show anything that is not significant 

# Trying to work with the response ratio data 
# I want to show by pre and post thaw (control and treatment) btt the res data does not specify that
ggplot(tax_response %>% slice(1:100), aes(x= ID, y = RR, color = Class)) + 
  geom_pointrange(aes(ymin = RR - err, ymax = RR + err)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw()
  






