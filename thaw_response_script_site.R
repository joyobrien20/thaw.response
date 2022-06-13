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
library(readxl)
library(ARPobservation) # checking the response ratio with another package
library(SingleCaseES)# checking the response ratio with another package pt. 2

# Installing and using mctools, it wasn't as straight forward as install.packages
# install.packages("devtools")
# devtools::install_github("leffj/mctoolsr")

# RESPONSE RATIO ANALYSIS
# Make a dataframe from the phyloseq object
dormOTU <- dorm1rarefied_OTU %>%
  data.frame()
# Read in file (load OTU table) 
thaw_response <- dormOTU
# We are going to scale each OTU count up by 0.01 so that 0 becomes 0.01 and the relatoionship stays the same
thaw_response <- thaw_response + 0.01

write.csv(thaw_response, "~/Desktop/thaw_response.csv")
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

# Response ratio calculations
# Sample variance
v1 = sd1^2/(ave1^2*num1); # Control #
v2 = sd2^2/(ave2^2*num2) # treatment #

# Square root of variance 
sqv = sqrt(v1+v2)

# RR = Response ratio
RR = log(ave2/ave1)
 
rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
  # I am going to move forward with RR (cited and vetted)

 
# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

ICl_R = RR-1.96*sqv
ICu_R = RR+1.96*sqv

ICl_r = rr-1.96*sqv
ICu_r = rr+1.96*sqv

# Error 
err_R = ICu_R-RR
err_r = ICu-rr

# Significance
Sig_R = ICl_R*ICu_R
Sig_adjust_R <- p.adjust(Sig_R, method = "fdr") # only 1 asv is significant ?? 

Sig_r = ICl_r*ICu_r
Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
R = as.data.frame(cbind(RR,err_R,Sig_R, ICl_R, ICu_R))
r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

R_adjust <- as.data.frame(cbind(RR,err_R,Sig_adjust_R, ICl_R, ICu_R))
r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))
# Create a new column from row names 
res$ID <- rownames(res)

# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# Come back to this and decide what to do when Hannah gets back to you
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

write.table(res,"thaw_response.txt",sep = "\t")

# Adding in the tax table information so that we can make sense of the res data frame
tax_response <- left_join(x = res, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(tax_response,"~/Desktop/tax_response.csv")
# Filter tax response so that we don't show anything that is not significant 

# Visualizing the RR data
ggplot(tax_response %>% slice(1:100), aes(x= ID, y = RR, color = Class)) + 
  geom_pointrange(aes(ymin = RR - err, ymax = RR + err)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw()
  


# RESPONSE RATIO ANALYSIS BY SITE! 
# CRREL

crrel_data <- extract(dormOTU, "X2_S73_L002", "X25_S96_L002", "X26_S97_L002",	"X27_S98_L002",	"X28_S99_L002",	"X29_S100_L002",	"X31_S102_L002",	"X4_S75_L002",	"X5_S76_L002",	"X6_S77_L002",	"X4_S75_L002",	"X5_S76_L002", "X6_S77_L002") 

crrel <- read_excel("~/Desktop/response_site.xlsx", sheet = "crrel")
# Make a dataframe from the phyloseq object
dormOTU <- dorm1rarefied_OTU %>%
 data.frame()

crrelOTU <- crrel %>%
  data.frame()
# Read in file (load OTU table) 
crrel_response <- as.numeric(crrelOTU)

# We are going to scale each OTU count up by 0.01 so that 0 becomes 0.01 and the relatoionship stays the same
crrelthaw_response <- crrel_response + 0.01


# Read in the CRREL OTU table that has been adjusted already (see above)

crrel <- read_excel("~/Desktop/response_site.xlsx", sheet = "crrel")

# Replace all 0 with NA
#thaw_response[thaw_response == 0] = 0.01
# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_crrel <- crrel;
ToF_crrel[ToF_crrel > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##

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

# Response ratio calculations
# Sample variance
v1 = sd1^2/(ave1^2*num1); # Control #
v2 = sd2^2/(ave2^2*num2) # treatment #

# Square root of variance 
sqv = sqrt(v1+v2)

# RR = Response ratio
RR = log(ave2/ave1)

rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
# I am going to move forward with RR (cited and vetted)


# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

ICl_R = RR-1.96*sqv
ICu_R = RR+1.96*sqv

ICl_r = rr-1.96*sqv
ICu_r = rr+1.96*sqv

# Error 
err_R = ICu_R-RR
err_r = ICu-rr

# Significance
Sig_R = ICl_R*ICu_R
Sig_adjust_R <- p.adjust(Sig_R, method = "fdr") # only 1 asv is significant ?? 

Sig_r = ICl_r*ICu_r
Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
R = as.data.frame(cbind(RR,err_R,Sig_R, ICl_R, ICu_R))
r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

R_adjust <- as.data.frame(cbind(RR,err_R,Sig_adjust_R, ICl_R, ICu_R))
r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))
# Create a new column from row names 
res$ID <- rownames(res)

# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# Come back to this and decide what to do when Hannah gets back to you
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

write.table(res,"thaw_response.txt",sep = "\t")

# Adding in the tax table information so that we can make sense of the res data frame
tax_response <- left_join(x = res, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(tax_response,"~/Desktop/tax_response.csv")
# Filter tax response so that we don't show anything that is not significant 

# Visualizing the RR data
ggplot(tax_response %>% slice(1:100), aes(x= ID, y = RR, color = Class)) + 
  geom_pointrange(aes(ymin = RR - err, ymax = RR + err)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw()

