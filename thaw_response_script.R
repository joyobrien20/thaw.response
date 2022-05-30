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
tax_table <- read_delim("~/Desktop/JMO_masters_2022/taxtable.txt")
#metadata <- read_excel("~/Desktop/incubation_16S_v4.xlsx", sheet = "metadata_final")
OTU <-read_delim("~/Desktop/JMO_masters_2022/otu_table.txt")
thaw.stats <- load_taxa_table("~/Desktop/JMO_masters_2022/taxtable.txt", "~/Desktop/JMO_masters_2022/otu_table.txt") # I cannot get this to work!
mapping_file <- thaw.stats$map_loaded 

# Response ratio analysis ## need to edit with your own specs and file names! 

# read in file (load ASV table) #
LP.ASV.response <- LP.stats$data_loaded
# Replace all 0 with NA #
LP.ASV.response[LP.ASV.response == 0] = NA
# creating a new varriable (TOF) new table w/ 0 = NA ###
ToF <- LP.ASV.response;
ToF[ToF > 0] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using realitive abundance table, creates average)
# treatment_names <- c("sample1", "sample2"...)
ave1 = apply(LP.ASV.response[,1:6],1,mean,na.rm = TRUE);# Control #
ave2 = apply(LP.ASV.response[,7:12],1,mean,na.rm = TRUE)# Treatment #
# Standard Deviation of ALL samples
sd1 = apply(LP.ASV.response[,1:6],1,sd,na.rm = TRUE);# Control #
sd2 = apply(LP.ASV.response[,7:12],1,sd,na.rm = TRUE)# treatment # 
# How many incounters for ASV/ number of observations ##
num1 = apply(ToF[,1:6],1,sum,na.rm = TRUE); #Control#
num2 = apply(ToF[,7:12],1,sum,na.rm = TRUE) # Treatment #
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
res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

write.table(res,"response_npkcm.txt",sep = "\t")
ggplot()