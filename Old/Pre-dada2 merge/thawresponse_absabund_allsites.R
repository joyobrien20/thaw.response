# Joy O'Brien
# Thaw response for all sites using absolute abundance data
library(mctoolsr)
library(ggplot2)
library(tibble)
library(dplyr)
library(multcompView)
library(phyloseq)
library(readxl)
library(biomformat)
library(readr)
library(readxl)
library(ARPobservation) # checking the response ratio with another package
library(SingleCaseES)# checking the response ratio with another package pt. 2
#Trying out [absolute abundance] here, so the paths have changed
absabund_RR <- absabund_allsites

# Getting rid of the column of numbers in the matrix <- crrel %>%
absabund_RR <- absabund_RR %>% 
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')

# We are going to scale each OTU count up by 0.01 so that 0 becomes 0.01 and the relationship stays the same
absabund_RR[absabund_RR == 0] = 0.01

# Save
#write.csv(absabund_RR, "~/Desktop/absabund_thawresponse_reseq_allsite.csv")


# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_absabund_RR <- absabund_RR;
ToF_absabund_RR[ToF_absabund_RR > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##


# Using relative abundance [absolute abundance] table, creates average
# treatment_names <- c("sample1", "sample2"...)
absabundave1 = apply(absabund_RR[,c(4,7:9,12,16:21,30:36)],1,mean,na.rm = TRUE);# Control 
absabundave2 = apply(absabund_RR[,c(1:3,5,6,10,11,13:15,22:29,37)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
absabundsd1 = apply(absabund_RR[,c(4,7:9,12,16:21,30:36)],1,sd,na.rm = TRUE);# Control #
absabundsd2 = apply(absabund_RR[,c(1:3,5,6,10,11,13:15,22:29,37)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
absabundnum1 = apply(ToF_absabund_RR[,c(4,7:9,12,16:21,30:36)],1,sum,na.rm = TRUE); #Control#
absabundnum2 = apply(ToF_absabund_RR[,c(1:3,5,6,10,11,13:15,22:29,37)],1,sum,na.rm = TRUE) # Treatment

# Response ratio calculations
# Sample variance
absabundv1 = absabundsd1^2/(absabundave1^2*absabundnum1); # Control #
absabundv2 = absabundsd2^2/(absabundave2^2*absabundnum2) # treatment #

# Square root of variance 
absabundsqv = sqrt(absabundv1+absabundv2)

# RR = Response ratio
absabundRR = log(absabundave2/absabundave1)

# Confidence interval

absabundICl_R = absabundRR-1.96*absabundsqv
absabundICu_R = absabundRR+1.96*absabundsqv


# Error 
absabunderr_R = absabundICu_R-absabundRR


# Significance
absabundSig_R = absabundICl_R*absabundICu_R
absabundSig_adjust_R <- p.adjust(absabundSig_R, method = "fdr") 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
allsites_R_abs = as.data.frame(cbind(absabundRR,absabunderr_R,absabundSig_R, absabundICl_R, absabundICu_R))

# Creating data frame of the adjusted values
allsites_R_adjust_abs <- as.data.frame(cbind(absabundRR,absabunderr_R,absabundSig_adjust_R, absabundICl_R, absabundICu_R))


# Create a new column from row names 
allsites_R_adjust_abs$ID <- rownames(allsites_R_adjust_abs)
#*****
# selecting rows (ASV's) with significance & values #
res = allsites_R_adjust_abs[absabundSig_adjust_R < 0.05,]

#write.table(allsites_R_adjust_abs,"crrel_abs_thaw_response.txt",sep = "\t")
allsites_R_adjust_abs$ID <- rownames(allsites_R_adjust_abs)
# Adding in the tax table information so that we can make sense of the res data frame

# ERROR HERE NEED TO FIX SO WE HAVE ALL DATA TOGETHER
# Load in the tax table from the R folder, rename to tax table instead of incubtation_16S reseq
allsites_tax_response_abs <- left_join(x = allsites_R_adjust_abs, y = incubation_16S_reseq, by = c("ID" = "ASV_ID"))

# Save to desktop
write_csv(allsites_tax_response_abs,"~/Desktop/allsites_absabund_thawresponse.csv")

# Filter tax response so that we don't show anything that is not significant (?)

crrel_final <- read_excel("~/Desktop/crrel_final_response.xls")
# write.csv(crrel_tax_response, "~/Desktop/crrel_tax_response.csv"

# Visualizing the RR data
ggplot(allsites_tax_response_abs %>% slice(1:100), aes(x= ID, y = absabundRR, color = Class)) + 
  geom_pointrange(aes(ymin = absabundRR - absabunderr_R, ymax = absabundRR + absabunderr_R)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("Log RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size = 12))

# END