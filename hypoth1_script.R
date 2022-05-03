# Joy O'Brien
# Master's Research 
# May 2, 2022

# Work to find evidence for OR against hypothesis 1: Permafrost soils with a high amount of dissolved 
# C and N pre that will have a higher cumulative respiration rate 

#******************************************************************
# May 3rd TO DO:
# RE-RUN THE GLM CODE WITH THE DOC_POST (LOG 10) BY SITE 
# RUN THE GLM CODE WITH THE DOC_PRE (LOG 10) BY SITE
# RUN THE GLM CODE WITH THE TDN_PRE (LOG10 IF NEEDED) BY SITE
# RUN THE GLM CODE WITH THE TDN_POST (LOG10 IF NEEDED) BY SITE
# ADD ALL TO THE PPT FOR LAB MEETING (MAY 6TH)
#*****************************************************************

#Load necessary library
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(e1071) #for skewness function

# Read in the data from excel
doctdn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "DOC_TDN")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum")

# Creating GLM (with the help of https://data-flair.training/blogs/generalized-linear-models-in-r/)

# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = doctdn$DOC_post) #change to DOC_pre, TDN_post, and TDN_pre

# We need to check if the dependent variable, DOC_pre in this case, is close to normal

par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post
plot(density(doctdn$DOC_post), main="DOC_post", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(doctdn$DOC_post), 3)))

polygon(density(doctdn$DOC_post), col="darkorange")

# Merge the data frames

tab<-merge(resp,doctdn)

rbind(resp, doctdn)
# creating the GLM #### DO THIS WITH EACH RESPONSE VARIABLE, AND NAME IT ACCORDINGLY, ALSO GROUP BY SITE-probably with the 
# actual running of the model and not in the code for ggplot 
LinearmodelA<-glm(Cumulative_Respiration ~ DOC_post, family = gaussian, data=tab)
print(LinearmodelA)
summary(LinearmodelA)

ggplot(LinearmodelA, aes(x= Cumulative_Respiration, y= DOC_post)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")

# creating LM
linearmodel1<- lm(Cumulative_Respiration ~ DOC_post, data=tab)
print(linearmodel1)
summary(linearmodel1)

ggplot(linearmodel1, aes(x= Cumulative_Respiration, y= DOC_post)) +
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")

model<- glm(Cumulative_Respiration ~ TDN, data = resp, doctdn, family = "binomial")


glimpse(doctdn)
glimpse(resp)

# Check for continuous variables
cont_doctdn <- select_if(doctdn, is.numeric)
cont_resp <- select_if(resp, is.numeric)
summary(cont_doctdn)
summary(cont_resp)

ggplot(cont_doctdn, aes(x = Sample_ID)) +
  geom_density(alpha = .2, fill = "#FF6666")

