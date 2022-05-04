# Joy O'Brien
# Master's Research 
# May 2, 2022

# Work to find evidence for OR against hypothesis 1: Permafrost soils with a high amount of dissolved 
# C and N pre that will have a higher cumulative respiration rate 

#Load necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(e1071) #for skewness function

# Read in the data from excel
doctdn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "DOC_TDN")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum")

# Creating GLM (with the help of https://data-flair.training/blogs/generalized-linear-models-in-r/)

# GLM for RESP and DOC_POST
# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = doctdn$DOC_post) #change to DOC_pre, TDN_post, and TDN_pre

# We need to check if the dependent variable, DOC_post in this case, is close to normal
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
respdoctdn<-merge(resp,doctdn)

# Creating the GLM #### DO THIS WITH EACH RESPONSE VARIABLE, AND NAME IT ACCORDINGLY, ALSO GROUP BY SITE-probably with the 
# actual running of the model and not in the code for ggplot ??
linearmodel_respDOCpost <- glm(Cumulative_Respiration ~ DOC_post, family = gaussian, data = respdoctdn)
# Look at the results of the glm
print(linearmodel_respDOCpost)
summary(linearmodel_respDOCpost)

# Let's see if we get an R2 value from running lm instead of glm
linearmodel_respDOCpost_test <- lm(Cumulative_Respiration ~ DOC_post, family = gaussian, data = respdoctdn)
print(linearmodel_respDOCpost_test)
summary(linearmodel_respDOCpost_test)
# Visualize the results of the glm with ggplot, note the log10 transformation
ggplot(linearmodel_respDOCpost, aes(x= Cumulative_Respiration, y = DOC_post) + #color = linearmodel_respDOCpost[["data"]][["site"]])) +
  geom_point() + # (aes(color = linearmodel_respDOCpost[["data"]][["site"]])) +
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm")

# This code below seems to work for getting one line through the sites 
ggplot(linearmodel_respDOCpost, aes(x = Cumulative_Respiration, y = DOC_post)) +
    geom_point(aes(color = linearmodel_respDOCpost[["data"]][["site"]])) +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm")

# Let's look at the plot that is not log transformed
ggplot(linearmodel_respDOCpost, aes(x = Cumulative_Respiration, y = DOC_post)) +
           geom_point() + 
           geom_smooth(method = "lm")

# EXTRA BELOW
# Check for continuous variables
#cont_doctdn <- select_if(doctdn, is.numeric)
#cont_resp <- select_if(resp, is.numeric)
#summary(cont_doctdn)
#summary(cont_resp)

#ggplot(cont_doctdn, aes(x = Sample_ID)) +
 # geom_density(alpha = .2, fill = "#FF6666")
#******************************************************
# GLM FOR RESP AND DOC PRE

# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = doctdn$DOC_pre) #change to DOC_pre, TDN_post, and TDN_pre

# We need to check if the dependent variable, DOC_pre in this case, is close to normal

par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post
plot(density(doctdn$DOC_post), main="DOC_pre", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(doctdn$DOC_pre), 3)))

polygon(density(doctdn$DOC_post), col="darkorange")

# creating the GLM #### DO THIS WITH EACH RESPONSE VARIABLE, AND NAME IT ACCORDINGLY, ALSO GROUP BY SITE-probably with the 
# actual running of the model and not in the code for ggplot 

linearmodel_respDOCpre <- glm(Cumulative_Respiration ~ DOC_pre, family = gaussian, data=respdoctdn)

print(linearmodel_respDOCpre)
summary(linearmodel_respDOCpre)

ggplot(linearmodel_respDOCpre, aes(x = Cumulative_Respiration, y = DOC_pre)) +
  geom_point(aes(color = linearmodel_respDOCpre[["data"]][["site"]])) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")

#**************************************************************

# GLM FOR RESP AND TDN POST

# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = doctdn$TDN_post) 

# We need to check if the dependent variable, DOC_pre in this case, is close to normal
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post
plot(density(doctdn$TDN_post), main="TDN_post", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(doctdn$TDN_post), 3)))

polygon(density(doctdn$TDN_post), col="darkorange")


# creating the GLM #### DO THIS WITH EACH RESPONSE VARIABLE, AND NAME IT ACCORDINGLY, ALSO GROUP BY SITE-probably with the 
# actual running of the model and not in the code for ggplot 
linearmodel_respTDNpost<-glm(Cumulative_Respiration ~ TDN_post, family = gaussian, data=respdoctdn)
print(linearmodel_respTDNpost)
summary(linearmodel_respTDNpost)

ggplot(linearmodel_respTDNpost, aes(x= Cumulative_Respiration, y= TDN_post)) +
  geom_point(aes(color=linearmodel_respTDNpost[["data"]][["site"]]))+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method="lm")

#**************************************************
#GLM FOR RESP AND TDN PRE

# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = doctdn$TDN_post) #change to DOC_pre, TDN_post, and TDN_pre

# We need to check if the dependent variable, DOC_pre in this case, is close to normal

par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of DOC_pre/post or TDN_pre/post
plot(density(doctdn$TDN_pre), main="TDN_pre", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(doctdn$DOC_post), 3)))

polygon(density(doctdn$TDN_pre), col="darkorange")

# creating the GLM #### DO THIS WITH EACH RESPONSE VARIABLE, AND NAME IT ACCORDINGLY, ALSO GROUP BY SITE-probably with the 
# actual running of the model and not in the code for ggplot 
linearmodel_respTDNpre<-glm(Cumulative_Respiration ~ TDN_pre, family = gaussian, data=respdoctdn)
print(linearmodel_respTDNpre)
summary(linearmodel_respTDNpre)

ggplot(linearmodel_respTDNpre, aes(x= Cumulative_Respiration, y= TDN_pre)) +
  geom_point(aes(color=linearmodel_respTDNpre[["data"]][["site"]]))+
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
