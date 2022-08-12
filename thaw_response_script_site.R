# Joy O'Brien
# May 30, 2022
# Master's project, Ernakovich Lab

# Thaw response code (response ratio analysis code from Hannah Holland-Moritz)
# This code is implemented to identify taxa  that respond positively (and negatively) to the thaw 
# disturbance. It is separated out as follows: all sites, CRREL, FL (Farmers Loop), and UT (Utqiagvik)

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

# RESPONSE RATIO ANALYSIS **** ALL SITES ****
# Make a dataframe from the phyloseq object
dormOTU <- dorm1rarefied_OTU %>%
  data.frame()

# Read in file (load OTU table) 
thaw_response <- dormOTU

# We are going to scale each OTU count up by 0.01 so that 0 becomes 0.01 and the relatoionship stays the same
thaw_response <- thaw_response + 0.01

# Save
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

# How many encounters for ASV number of observations
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
 
#rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
  # I am going to move forward with RR (cited and vetted)

 
# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval
ICl_R = RR-1.96*sqv
ICu_R = RR+1.96*sqv

#ICl_r = rr-1.96*sqv
#ICu_r = rr+1.96*sqv

# Error 
err_R = ICu_R-RR
#err_r = ICu-rr

# Significance
Sig_R = ICl_R*ICu_R
Sig_adjust_R <- p.adjust(Sig_R, method = "fdr") # adjusting for multiple comparisons

#Sig_r = ICl_r*ICu_r
#Sig_adjust_r <- p.adjust(Sig_r, method = "fdr")

# Creating data frame w/ 3 coloums (RR, err, Sig) 
R = as.data.frame(cbind(RR,err_R,Sig_R, ICl_R, ICu_R))
#r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

R_adjust <- as.data.frame(cbind(RR,err_R,Sig_adjust_R, ICl_R, ICu_R))
#r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))
# Create a new column from row names 
R_adjust$ID <- rownames(R_adjust)

# selecting rows (ASV's) with significance & values #
# res = res[res$Sig > 0,]
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

write.table(R_adjust,"thaw_response.txt",sep = "\t")

# Adding in the tax table information so that we can make sense of the res data frame
tax_response <- left_join(x = R_adjust, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(tax_response,"~/Desktop/tax_response.csv")
# Filter tax response so that we don't show anything that is not significant 

# Visualizing the RR data
ggplot(tax_response %>% slice(1:100), aes(x= ID, y = RR, color = Class)) + 
  geom_pointrange(aes(ymin = RR - err, ymax = RR + err)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw()

# Final visualization
ggplot(tax_response %>% slice(1:100), aes(x= ID, y = RR, color = Class)) + 
  geom_pointrange(aes(ymin = RR - err_R, ymax = RR + err_R)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  theme(text = element_text(size = 12))
  
#****************************************************
# RESPONSE RATIO ANALYSIS BY SITE! 
# CRREL

#Trying out [absolute abundance] here, so the paths have changed
crrel <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "CRREL_final_absabund")

# Getting rid of the column of numbers in the matrix
crrel <- crrel %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')

# We are going to scale each OTU count up by 0.01 so that 0 becomes 0.01 and the relationship stays the same
crrel[crrel == 0] = 0.01

# Save
write.csv(crrel, "~/Desktop/CRREL_abs_thawresponse.csv")


# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_crrel <- crrel;
ToF_crrel[ToF_crrel > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using relative abundance [absolute abundance] table, creates average
# treatment_names <- c("sample1", "sample2"...)
crrelave1 = apply(crrel[,c(8:14)],1,mean,na.rm = TRUE);# Control 
crrelave2 = apply(crrel[,c(1:7)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
crrelsd1 = apply(crrel[,c(8:14)],1,sd,na.rm = TRUE);# Control #
crrelsd2 = apply(crrel[,c(1:7)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
crrelnum1 = apply(ToF_crrel[,c(8:14)],1,sum,na.rm = TRUE); #Control#
crrelnum2 = apply(ToF_crrel[,c(1:7)],1,sum,na.rm = TRUE) # Treatment #

# Response ratio calculations
# Sample variance
crrelv1 = crrelsd1^2/(crrelave1^2*crrelnum1); # Control #
crrelv2 = crrelsd2^2/(crrelave2^2*crrelnum2) # treatment #

# Square root of variance 
crrelsqv = sqrt(crrelv1+crrelv2)

# RR = Response ratio
crrelRR = log(crrelave2/crrelave1)

#rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
# I am going to move forward with RR (cited and vetted)


# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

crrelICl_R = crrelRR-1.96*crrelsqv
crrelICu_R = crrelRR+1.96*crrelsqv

#ICl_r = rr-1.96*sqv
#ICu_r = rr+1.96*sqv

# Error 
crrelerr_R = crrelICu_R-crrelRR
#err_r = ICu-rr

# Significance
crrelSig_R = crrelICl_R*crrelICu_R
crrelSig_adjust_R <- p.adjust(crrelSig_R, method = "fdr") 

#Sig_r = ICl_r*ICu_r
#Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
crrel_R_abs = as.data.frame(cbind(crrelRR,crrelerr_R,crrelSig_R, crrelICl_R, crrelICu_R))
#r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

crrel_R_adjust_abs <- as.data.frame(cbind(crrelRR,crrelerr_R,crrelSig_adjust_R, crrelICl_R, crrelICu_R))
#r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))

# Create a new column from row names 
crrel_R_adjust_abs$ID <- rownames(crrel_R_adjust_abs)

# selecting rows (ASV's) with significance & values #
crrel_R_adjust_abs = crrel_R_adjust_abs[crrel_R_adjust_abs$Sig > 0,]
# Come back to this and decide what to do when Hannah gets back to you
# crrel_R_adjust_abs_filter <- filter(crrel_R_adjust_abs, Sig > 0)

# creating new coloumn named "ID" # 
crrel_R_adjust_abs[,4] = rownames(crrel_R_adjust_abs);
names(crrel_R_adjust_abs)[4] = "ID"

write.table(crrel_R_adjust_abs,"crrel_abs_thaw_response.txt",sep = "\t")
crrel_R_adjust_abs$ID <- rownames(crrel_R_adjust_abs)
# Adding in the tax table information so that we can make sense of the res data frame
crrel_tax_response_abs <- left_join(x = crrel_R_adjust_abs, y = tax_tablecsv, by = c("ID" = "#ASV_ID"))

# Save to desktop
write_csv(crrel_tax_response_abs,"~/Desktop/crrel)tax_response_final.csv")

# Filter tax response so that we don't show anything that is not significant (?)

crrel_final <- read_excel("~/Desktop/crrel_final_response.xls")
# write.csv(crrel_tax_response, "~/Desktop/crrel_tax_response.csv"

# Visualizing the RR data
ggplot(crrel_tax_response_abs %>% slice(1:100), aes(x= ID, y = crrelRR, color = Class)) + 
  geom_pointrange(aes(ymin = crrelRR - crrelerr_R, ymax = crrelRR + crrelerr_R)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("Log RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  theme(text = element_text(size = 12))
  
#*********************************************************************************
 
# FARMERS LOOP

# For absolute abundance
FL <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "FL_final_absabund")
#FL <- read_excel("~/Desktop/response_site.xlsx", sheet = "FL")

FL[FL == 0] = 0.01
# Getting rid of the column of numbers in the matrix
FL <- FL %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')

# Replace all 0 with NA
#thaw_response[thaw_response == 0] = 0.01
# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_FL <- FL;
ToF_FL[ToF_FL > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using relative abundance table, creates average
# treatment_names <- c("sample1", "sample2"...)
FLave1 = apply(FL[,c(8:14)],1,mean,na.rm = TRUE);# Control 
FLave2 = apply(FL[,c(1:7)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
FLsd1 = apply(FL[,c(8:14)],1,sd,na.rm = TRUE);# Control #
FLsd2 = apply(FL[,c(1:7)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
FLnum1 = apply(ToF_FL[,c(8:14)],1,sum,na.rm = TRUE); #Control#
FLnum2 = apply(ToF_FL[,c(1:7)],1,sum,na.rm = TRUE) # Treatment #

# Response ratio calculations
# Sample variance
FLv1 = FLsd1^2/(FLave1^2*FLnum1); # Control #
FLv2 = FLsd2^2/(FLave2^2*FLnum2) # treatment #

# Square root of variance 
FLsqv = sqrt(FLv1+FLv2)

# RR = Response ratio
FLRR = log(FLave2/FLave1)

#rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
# I am going to move forward with RR (cited and vetted)


# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

FLICl_R = FLRR-1.96*FLsqv
FLICu_R = FLRR+1.96*FLsqv

#ICl_r = rr-1.96*sqv
#ICu_r = rr+1.96*sqv

# Error 
FLerr_R = FLICu_R-FLRR
#err_r = ICu-rr

# Significance
FLSig_R = FLICl_R*FLICu_R
FLSig_adjust_R <- p.adjust(FLSig_R, method = "fdr") 

#Sig_r = ICl_r*ICu_r
#Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
FL_R = as.data.frame(cbind(FLRR,FLerr_R,FLSig_R, FLICl_R, FLICu_R))
#r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

FL_R_adjust <- as.data.frame(cbind(FLRR,FLerr_R,FLSig_adjust_R, FLICl_R, FLICu_R))
#r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))

# Create a new column from row names 
FL_R_adjust$ID <- rownames(FL_R_adjust)

# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

#write.table(res,"thaw_response.txt",sep = "\t")
# Adding in the tax table information so that we can make sense of the res data frame
FL_abs_tax_response <- left_join(x = FL_R_adjust, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(FL_abs_tax_response,"~/Desktop/FL_abs_tax_response_final.csv")


# Visualizing the RR data
ggplot(FL_abs_tax_response %>% slice(1:100), aes(x= ID, y = FLRR, color = Class)) + 
  geom_pointrange(aes(ymin = FLRR - FLerr_R, ymax = FLRR + FLerr_R)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("Log RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  theme(text = element_text(size = 12))

#******************************************************
#*Utqiagvik response ratio 

# For absolute abundance

UT <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "UT_final_absabund")
# Relative abundance 
# UT <- read_excel("~/Desktop/response_site.xlsx", sheet = "UT")

# Getting rid of the column of numbers in the matrix
UT <- UT %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')

UT[UT == 0] = 0.01
# Replace all 0 with NA
#thaw_response[thaw_response == 0] = 0.01
# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_UT <- UT;
ToF_UT[ToF_UT > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using relative abundance table, creates average
# treatment_names <- c("sample1", "sample2"...) 
UTave1 = apply(UT[,c(9:16)],1,mean,na.rm = TRUE);# Control
UTave2 = apply(UT[,c(1:8)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
UTsd1 = apply(UT[,c(1:8)],1,sd,na.rm = TRUE);# Control #
UTsd2 = apply(UT[,c(9:16)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
UTnum1 = apply(ToF_UT[,c(1:8)],1,sum,na.rm = TRUE); #Control#
UTnum2 = apply(ToF_UT[,c(9:16)],1,sum,na.rm = TRUE) # Treatment #

# Response ratio calculations
# Sample variance
UTv1 = UTsd1^2/(UTave1^2*UTnum1); # Control #
UTv2 = UTsd2^2/(UTave2^2*UTnum2) # treatment #

# Square root of variance 
UTsqv = sqrt(UTv1+UTv2)

# RR = Response ratio
UTRR = log(UTave2/UTave1)

#rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
# I am going to move forward with RR (cited and vetted)


# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

UTICl_R = UTRR-1.96*UTsqv
UTICu_R = UTRR+1.96*UTsqv

#ICl_r = rr-1.96*sqv
#ICu_r = rr+1.96*sqv

# Error 
UTerr_R = UTICu_R-UTRR
#err_r = ICu-rr

# Significance
UTSig_R = UTICl_R*UTICu_R
UTSig_adjust_R <- p.adjust(UTSig_R, method = "fdr") 

#Sig_r = ICl_r*ICu_r
#Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
UT_R = as.data.frame(cbind(UTRR,UTerr_R,UTSig_R, UTICl_R, UTICu_R))
#r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

UT_R_adjust <- as.data.frame(cbind(UTRR,UTerr_R,UTSig_adjust_R, UTICl_R, UTICu_R))
#r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))

# Create a new column from row names 
UT_R_adjust$ID <- rownames(UT_R_adjust)

# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

#write.table(res,"thaw_response.txt",sep = "\t")
# Adding in the tax table information so that we can make sense of the res data frame
UT_abs_tax_response <- left_join(x = UT_R_adjust, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(UT_abs_tax_response,"~/Desktop/UT_abs_tax_response_final.csv")

# Visualizing the RR data
ggplot(UT_abs_tax_response %>% slice(1:100), aes(x= ID, y = UTRR, color = Class)) + 
  geom_pointrange(aes(ymin = -40, ymax = 40)) + # had to change this because the scale was too big and none of the points were visible
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("Log RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  theme(text = element_text(size = 12))

# ABSOLUTE ABUNDANCE FOR ALL SITES

ALL <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "abs_abund_all")
# Relative abundance 
# UT <- read_excel("~/Desktop/response_site.xlsx", sheet = "UT")

# Getting rid of the column of numbers in the matrix
ALL <- ALL %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ASV')

ALL[ALL == 0] = 0.01
# Replace all 0 with NA
#thaw_response[thaw_response == 0] = 0.01
# creating a new variable (TOF) new table w/ 0 = NA ###
ToF_ALL <- ALL;
ToF_ALL[ToF_ALL > 0.01] <- 1 #anything above 0 is 1 (Presence absence table) ##

# Using relative abundance table, creates average
# treatment_names <- c("sample1", "sample2"...) 
ALLave1 = apply(ALL[,c(8:14, 22:28, 37:44)],1,mean,na.rm = TRUE);# Control
ALLave2 = apply(ALL[,c(1:7, 15:21, 29:36)],1,mean,na.rm = TRUE)# Treatment 

# Standard Deviation of ALL samples
ALLsd1 = apply(ALL[,c(8:14, 22:28, 37:44)],1,sd,na.rm = TRUE);# Control #
ALLsd2 = apply(ALL[,c(1:7, 15:21, 29:36)],1,sd,na.rm = TRUE)# treatment # 

# How many encounters for ASV/ number of observations
ALLnum1 = apply(ToF_ALL[,c(8:14, 22:28, 37:44)],1,sum,na.rm = TRUE); #Control#
ALLnum2 = apply(ToF_ALL[,c(1:7, 15:21, 29:36)],1,sum,na.rm = TRUE) # Treatment #

# Response ratio calculations
# Sample variance
ALLv1 = ALLsd1^2/(ALLave1^2*ALLnum1); # Control #
ALLv2 = ALLsd2^2/(ALLave2^2*ALLnum2) # treatment #

# Square root of variance 
ALLsqv = sqrt(ALLv1+ALLv2)

# RR = Response ratio
ALLRR = log(ALLave2/ALLave1)

#rr <- log(ave2 + ave1)-log(ave2) #this was in an effort to compare response ratio calculations +
# I am going to move forward with RR (cited and vetted)


# logRespRatio(aves, phase = sqv, base_level = 0.01)

# Confidence interval

ALLICl_R = ALLRR-1.96*ALLsqv
ALLICu_R = ALLRR+1.96*ALLsqv

#ICl_r = rr-1.96*sqv
#ICu_r = rr+1.96*sqv

# Error 
ALLerr_R = ALLICu_R-ALLRR
#err_r = ICu-rr

# Significance
ALLSig_R = ALLICl_R*ALLICu_R
ALLSig_adjust_R <- p.adjust(ALLSig_R, method = "fdr") 

#Sig_r = ICl_r*ICu_r
#Sig_adjust_r <- p.adjust(Sig_r, method = "fdr") # only 1 asv is significant ?? 

# Creating data frame w/ 3 coloums (RR, err, Sig) 
ALL_R = as.data.frame(cbind(ALLRR,ALLerr_R,ALLSig_R, ALLICl_R, ALLICu_R))
#r = as.data.frame(cbind(rr,err_r,Sig_r, ICl_r, ICu_r))

ALL_R_adjust <- as.data.frame(cbind(ALLRR,ALLerr_R,ALLSig_adjust_R, ALLICl_R, ALLICu_R))
#r_a = as.data.frame(cbind(rr,err_r,Sig_adjust_r, ICl_r, ICu_r))

# Create a new column from row names 
ALL_R_adjust$ID <- rownames(ALL_R_adjust)

# selecting rows (ASV's) with significance & values #
#res = res[res$Sig > 0,]
# res <- filter(res, Sig > 0)

# creating new coloumn named "ID" # 
#res[,4] = rownames(res);
#names(res)[4] = "ID"

#write.table(res,"thaw_response.txt",sep = "\t")
# Adding in the tax table information so that we can make sense of the res data frame
ALL_abs_tax_response <- left_join(x = ALL_R_adjust, y = tax_tablecsv, by = c("ID"="#ASV_ID"))

# Save to desktop
write_csv(ALL_abs_tax_response,"~/Desktop/ALL_abs_tax_response_final.csv")

# Visualizing the RR data
ggplot(ALL_abs_tax_response %>% slice(1:100), aes(x= ID, y = ALLRR, color = Class)) + 
  geom_pointrange(aes(ymin = ALLRR - ALLerr_R, ymax = ALLRR + ALLerr_R)) + # had to change this because the scale was too big and none of the points were visible
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~Class,scales = "free_x") +
  theme_bw() +
  ylab("Log RR") +
  xlab ("ASV ID") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  theme(text = element_text(size = 12))
# END 
