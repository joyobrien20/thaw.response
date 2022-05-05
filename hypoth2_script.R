# Joy O'Brien
# Master's research
# May 3, 2022

# Work to find evidence for or against hypothesis 2: permafrost soils with a higher cumulative respiration rate 
#will have a higher qPCR abundance post thaw
#***********************************************
# NOTES FOR MAY 5TH
# CLEAN UP LAB MEETING SLIDES 
# LOOK MORE INTO THE OUTPUT RESULTS OF GLM AND WHAT THEY MEAN, ALSO LOOK INTO DIFFERENT TYPES OF DISTRIBUTION FUNCTIONS 
# RUN YOUR ANOVA COMPARISONS IN THE ANOVA SCRIPT (may have to alter/set up the excel files for this ) 
# For hypoth 3: create a new script, 
#***********************************
# ##### NOTE TO SELF: COME BACK UP HERE AND CHECK ALL SITES WITH THE X = COPY NUMBER, Y = RESPIRATION !
# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(e1071) #for skewness function

# X = COPY NUMBERS, Y= CUMULATIVE RESPIRATION: CODE BELOW******************************************
#NOTE: EXTREME OUTLIER, AT4 REP 8 WAS REMOVED FROM THE FOLLOWING DATA SET 
# Read in the data from excel
cpnumb <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "copy_numbers")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum (2)")

# Visualizing the data to start off 
scatter.smooth(cpnumb$copy_number_per_gsoil_PRE, y = resp$Cumulative_Respirationcpnumb)

# Checking for normality--totally forgot to do this! 
shapiro.test(cpnumb$copy_number_per_gsoil_PRE)
shapiro.test(resp$Cumulative_Respiration)

# We need to check if the dependent variable, DOC_post in this case, is close to normal
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of copy numbers pre
plot(density(cpnumb$copy_number_per_gsoil_PRE), main="copy_number_per_gsoil_PRE", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(cpnumb$copy_number_per_gsoil_PRE), 3)))

polygon(density(cpnumb$copy_number_per_gsoil_PRE), col="darkorange")

# Merge the data frames
respcpnumb <- merge(resp,cpnumb) 

# Creating the glm
linearmodel_respcopypre <- glm(copy_number_per_gsoil_PRE ~ Cumulative_Respiration, family=gaussian, data = respcpnumb)

# Look at the results of the glm
print(linearmodel_respcopypre)
summary(linearmodel_respcopypre)

# Visualize with ggplot  ## Okay this is fine, BUT I want to show this by site so theoretically make 3 plots for pre and 3 for post 
ggplot(linearmodel_respcopypre, aes(x = copy_number_per_gsoil_PRE, y = Cumulative_Respiration)) +
  geom_point (aes(color = linearmodel_respcopypre[["data"]][["site"]])) + 
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm")

#**************************************************************************

# ANALYZING POST-THAW COPY NUMBERS AND RESPIRATION, X=COPY NUMBERS AND Y=RESPIRATION

# Visualizing the data to start off 
scatter.smooth(x = cpnumb$copy_number_per_gsoil_POST, y = resp$Cumulative_Respiration)

# Checking for normality in the post thaw copy number values 
shapiro.test(cpnumb$copy_number_per_gsoil_POST)

# We need to check if the dependent variable is close to normal (?)
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of copy numbers pre
plot(density(cpnumb$copy_number_per_gsoil_POST), main="copy_number_per_gsoil_POST", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(cpnumb$copy_number_per_gsoil_POST), 3)))

polygon(density(cpnumb$copy_number_per_gsoil_PRE), col="darkorange")

# Merge the data frames
respcpnumb <- merge(resp,cpnumb)

linearmodel_respcopypost <- glm(copy_number_per_gsoil_POST ~ Cumulative_Respiration, family=gaussian, data = respcpnumb)

# Look at the results of the glm
print(linearmodel_respcopypost)
summary(linearmodel_respcopypost)

# Visualize with ggplot # HAVE TO GET THIS BY SITE AGAIN, IT'S FREAKING OUT OTHERWISE 
ggplot(linearmodel_respcopypost, aes(x = copy_number_per_gsoil_POST, y = Cumulative_Respiration)) +
  geom_point(aes(color = linearmodel_respcopypost[["data"]][["site"]])) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")

#**********************************************************************
#BELOW IS THE SAME CODE AS ABOVE BUT WITH X=RESPIRATION, Y=COPY NUMBERS
#**********************************************************************

#*HERE IS WHERE YOU WANT TO LOOK AT IT BY SITE CAUSE IT LOOKS FLAT AND OFF IN THE PRE-LIM DATA CHECK AND THE GLM 

# Visualizing the data to start off 
scatter.smooth(x = resp$Cumulative_Respiration, y = cpnumb$copy_number_per_gsoil_PRE)

# We need to check if the dependent variable, DOC_post in this case, is close to normal
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of copy numbers pre
plot(density(cpnumb$copy_number_per_gsoil_PRE), main="copy_number_per_gsoil_PRE", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(cpnumb$copy_number_per_gsoil_PRE), 3)))

polygon(density(cpnumb$copy_number_per_gsoil_PRE), col="darkorange")

# Merge the data frames
respcpnumb <- merge(resp,cpnumb) 

# Creating the glm
linearmodel_respcopypre <- glm(copy_number_per_gsoil_PRE ~ Cumulative_Respiration, family=gaussian, data = respcpnumb)

# Look at the results of the glm
print(linearmodel_respcopypre)
summary(linearmodel_respcopypre)

# Visualize with ggplot  ## Okay this is fine, BUT I want to show this by site so theoretically make 3 plots for pre and 3 for post 
ggplot(linearmodel_respcopypre, aes(x = copy_number_per_gsoil_PRE, y = Cumulative_Respiration)) +
  geom_point (aes(color = linearmodel_respcopypre[["data"]][["site"]])) + 
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm")

#**************************************************************************

# ANALYZING POST-THAW COPY NUMBERS AND RESPIRATION, X=COPY NUMBERS AND Y=RESPIRATION

# Visualizing the data to start off 
scatter.smooth(x = cpnumb$copy_number_per_gsoil_POST, y = resp$Cumulative_Respiration)

# We need to check if the dependent variable is close to normal (?)
par(mfrow=c(1, 2)) # divide graph area in 2 columns

# Look at the density of cumulative respiration
plot(density(resp$Cumulative_Respiration), main="Cumulative_Respiration", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(resp$Cumulative_Respiration), 3)))
polygon(density(resp$Cumulative_Respiration), col="darkorange")

# Look at the density of copy numbers pre
plot(density(cpnumb$copy_number_per_gsoil_POST), main="copy_number_per_gsoil_POST", ylab="Frequency", 
     sub=paste("Skewness:", round(e1071::skewness(cpnumb$copy_number_per_gsoil_POST), 3)))

polygon(density(cpnumb$copy_number_per_gsoil_PRE), col="darkorange")

# Merge the data frames
respcpnumb <- merge(resp,cpnumb)

linearmodel_respcopypost <- glm(Cumulative_Respiration ~ copy_number_per_gsoil_POST, family=gaussian, data = respcpnumb)

# Look at the results of the glm
print(linearmodel_respcopypost)
summary(linearmodel_respcopypost)

# Visualize with ggplot # HAVE TO GET THIS BY SITE AGAIN, IT'S FREAKING OUT OTHERWISE, oKAY BUT IT WAS FREAKING OUT BECASUE OF THE OUTLIER REP IT DIDNT KNOW WHAT TO DO WITH IT 
ggplot(linearmodel_respcopypost, aes(x = Cumulative_Respiration, y = copy_number_per_gsoil_POST)) +
  geom_point(aes(color = linearmodel_respcopypost[["data"]][["site"]])) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")
