# Joy O'Brien
# Master's research
# May 3, 2022

# Work to find evidence for or against hypothesis 2: permafrost soils with a higher cumulative respiration rate 
#will have a higher qPCR abundance post thaw
#***********************************************
# NOTES FOR MAY 4TH
# RUN THE GLM ON POST THAW COPY NUMBERS AND RESPIRATION, INCLUDE THIS IN THE PPT
# LOOK MORE INTO THE OUTPUT RESULTS OF GLM AND WHAT THEY MEAN, ALSO LOOK INTO DIFFERENT TYPES OF DISTRIBUTION FUNCTIONS 
# START YOUR ANOVA COMPARISONS 
#***********************************
# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(e1071) #for skewness function

# Read in the data from excel
cpnumb <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "copy_numbers")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum (2)")

# Visualizing the data to start off 
scatter.smooth(x=resp$Cumulative_Respiration, y = cpnumb$copy_number_per_gsoil_PRE)

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
respcpnumb <- merge(resp,cpnumb) #This is not merging properly because I am ending up with 0 obs in the data frame

linearmodel_respcopypre <- glm(Cumulative_Respiration ~ copy_number_per_gsoil_PRE, family=gaussian, data = respcpnumb)

# Look at the results of the glm
print(linearmodel_respcopypre)
summary(linearmodel_respcopypre)

# Visualize with ggplot 
ggplot(linearmodel_respcopypre, aes(x = Cumulative_Respiration, y = copy_number_per_gsoil_PRE)) +
  geom_point(aes(color = linearmodel_respcopypre[["data"]][["site"]])) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")

