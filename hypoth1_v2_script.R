# Joy O'Brien
# Master's work 
# May 6, 2022
# This script is a re-work of hypoth1_script taking into account advice from lab meeting on May 6, 2022

# In this script we will work to: try a Spearman correlation b/w DOC and TDN and cumulative respiration
# Also we will be working forward with DOC and TDN on the x-axis and respiration on the y-axis

# TO DO'S:
#     Note: If you do a glm you need to possibly transform the data prior, and add the pre-values in and account for site 
#   Run a Spearman correlation with the difference in DOC and TDN
#   Run a GLM with the difference in DOC and the difference in TDN, with cumulative respiration as the y axis (response)
#   Confirm with Jessica that GLM or Spearman is the way to go; once confirmed and if needed, 
#   run glm of all comparisons with the family function not gaussian, and an added link function

#   If you are going to run a linear mixed effect model, make sure that site is the mixed effect and that others are the fixed effects
#   RUN SPEARMAN CORRELATION ON SITES WITH THE RESP, THE COPY NUMBERS, THE DOC AND THE TDN?

#Load necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(e1071) #for skewness function
library(qqplotr)

# Read in the data from excel
doctdn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "DOC_TDN_post")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum")
copynumb <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "copy_numbers")
resp4copy <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum_forcopynumb")

# Performing Spearman rank correlations on DOC_post vs. respiration

# Make a QQ plot 
qqnorm(x = doctdn$DOC_post, y = resp$Cumulative_Respiration, main = "Q-Q Plot")

# Trying a correlation test for DOC_post and respiration  
doccorr <- cor.test(x = doctdn$DOC_post, y = resp$Cumulative_Respiration, method = "spearman") #p value = 2.045e-05, rho = 0.766087
print(doccorr)

# Merge the data frames
respdoctdn <- merge(resp,doctdn)

# Visualize with ggplot
ggplot(respdoctdn, aes(x = DOC_post, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, shape = core)) +
  scale_x_log10() +
  scale_y_log10()

# Performing Spearman rank correlations on TDN_post and respiration

tdncorr <- cor.test(x = doctdn$TDN_post, y = resp$Cumulative_Respiration, method = "spearman")
print(tdncorr) # p-value = 0.001437, rho = 0.6234783

# Visualize with ggplot 
ggplot(respdoctdn, aes(x = TDN_post, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, shape = core)) +
  scale_x_log10() +
  scale_y_log10()

# This is the code that Hannah wrote to have the TDN_pre come in as dotted lines for comparison ##IF YOU ARE GOING TO USE THIS YOU NEED TO WORK THE PRE VALUES INTO THE FRAME YOU HAVE NOW 
ggplot(respdoctdn, aes(x = TDN_post, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, shape = core)) +
  geom_vline(aes(xintercept = TDN_pre, color = site), linetype = "dashed") +
  #geom_point(data = respdoctdn, aes(x = TDN_pre), color = "red") +
  stat_smooth(method = "lm") +
  scale_x_log10()

# LINEAR MIXED MODELING CODE BELOW
library(lme4)
respdoctdn$TDN_post.log <- log(respdoctdn$TDN_post)
tdn_resp.lme <- lme4::lmer(Cumulative_Respiration ~ TDN_post + 1|site), data = respdoctdn)

save("~/Desktop/joys_rdata.Rdata")
  
#************************************************************************
# PERFORMING A GLM ON DOC_POST VS CUMULATIVE RESPIRATION 

# Visualizing the data to start off 
scatter.smooth(x = doctdn$DOC_post, y = resp$Cumulative_Respiration) 

# We need to check if the dependent variable, DOC_post in this case, is close to normal
par(mfrow = c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main = "Cumulative_Respiration", ylab = "Frequency", 
     sub = paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col = "darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post
plot(density(doctdn$DOC_post), main = "DOC_post", ylab = "Frequency", 
     sub = paste("Skewness:", round(e1071::skewness(doctdn$DOC_post), 3)))

polygon(density(doctdn$DOC_post), col = "darkorange")

# Merge the data frames
respdoctdn <- merge(resp,doctdn)

# Visualizing with ggplot 
ggplot(respdoctdn, aes(x = DOC_post, y = Cumulative_Respiration)) + 
  geom_point()
  
# Creating the GLM
linearmodel_respDOCpost <- glm(DOC_post ~ Cumulative_Respiration, family = gaussian, data = respdoctdn) #this needs to be fixed and site needs to be taken into account in addition to the pre-thaw values

# Look at the results of the glm
print(linearmodel_respDOCpost)
summary(linearmodel_respDOCpost)

# Visualize the results of the glm with ggplot, note the log10 transformatioN
# This code below seems to work for getting one line through the sites
ggplot(linearmodel_respDOCpost, aes(x = DOC_post, y = Cumulative_Respiration)) +
  geom_point(aes(color = linearmodel_respDOCpost[["data"]][["site"]])) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")

#********************************************************************************
# GLM FOR RESP AND TDN POST (need to account for TDN pre within the model and site)

# Visualizing the data to start off 
scatter.smooth( x = doctdn$TDN_post, y=resp$Cumulative_Respiration) 

# We need to check if the dependent variable,TDN_ post is close to normal
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post (whichever section of the code you are in)
plot(density(doctdn$TDN_post), main="TDN_post", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(doctdn$TDN_post), 3)))

polygon(density(doctdn$TDN_post), col="darkorange")

# checking for normality 
shapiro.test(doctdn$TDN_post) #p-value = 0.008402, this is also not normally distributed

# creating the GLM # need to come back and make changes to this! 

linearmodel_respTDNpost<-glm(TDN_post ~ Cumulative_Respiration, family = gaussian, data=respdoctdn)
print(linearmodel_respTDNpost)
summary(linearmodel_respTDNpost)

ggplot(linearmodel_respTDNpost, aes(x = TDN_post, y = Cumulative_Respiration)) +
  geom_point(aes(color=linearmodel_respTDNpost[["data"]][["site"]]))+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")
     

     
     
                                    
                                   