# Joy O'Brien
# Master's work
# May 4, 2022

# THIS IS THE ANOVAS FOR SITE COMPARISON SCIPT
# 
# Data we are working with here: We will be running ANOVAs on DOC data post thaw, TDN data post thaw, copy numbers pre thaw and copy numbers post thaw and cum respiration across sites

# The reason we are NOT including pre thaw data here is because the pre thaw data is not properly replicated (due to lack of soil to replicate). 
# RUN A TUKEY TEST AFTER THE ANOVA TEST, NON PARAMETRIC ALTERNATIVE: SHAPIRO, KRUSKALL, DUNN 
# Additionally, % C, % N, C:N ratio, pH (average), TOC, and TON pre and post WILL BE IN A TABLE FORMAT (see notebook number 6 page 22 for reference)

# TO DO:
# NEED TO CHECK TO SEE IF YOU ARE RUNNING THE KRUSKALL WITH THE RIGHT PARAMETERS


# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(FSA) #installed for the dunntest

# Read in the data 
doctdn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "DOC_TDN")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum (2)")

# FIRST: RUNNING ANOVA ON DOC POST THAW

# Checking for normality with a Shapiro test 
shapiro.test(doctdn$DOC_post) #p-value = 0.02452, not normal

# The data is NOT normally distrubuted so I am going to try a Kruskall-Wallis test 
kruskal.test(DOC_post ~ site, data = doctdn) # Kruskal-Wallis chi-squared = 14.945, df = 2, p-value = 0.0005685

# Since the p-value of the kruskal test is significant, we are going to run a DUNN test 
# which I think is comparable to a tukey after an ANOVA 

dunnTest(DOC_post ~ site, data = doctdn, method = "bonferroni") # Check to see if the method is necessary? If it needs to be changed 

#let's visualize DOC_post content across site 
boxplot(DOC_post ~ site, data = doctdn)


# Visualizing/checking for theoretical model distribution (https://www.statology.org/anova-assumptions/)
#hist(doctdn$DOC_post)
#qqnorm(anova_docpost$residuals)

# NOW WE WILL BE MOVING ON TO ASSESS TND_POST BY SITE 
shapiro.test(doctdn$TDN_post) #p-value = 0.008402; IT IS NOT NORMALLY DISTRIBUTED 

# Therefore, we will run a Kruskall-Wallis test 
kruskal.test(TDN_post ~ site, data = doctdn) #chi-squared = 19.86, df = 2, p-value = 4.869e-05

# Since the p-value of the Kruskal for TDN_post is significant, we need to run a dunntest

dunnTest(TDN_post ~ site, data = doctdn, method = "bonferroni") 

# Now let's visualize with a boxplot
boxplot(TDN_post ~ site, data = doctdn)

# Now let's look at a dunn test
dunnTest(TDN_post ~ site, data = doctdn, method = "bonferroni")

# NOW TRANSITIONING TO CUMULATIVE RESPIRATION

# Let's look at the differences in cumulative respiration across sites 
shapiro.test(resp$Cumulative_Respiration) #this is not normally distributed W = 0.6816, p-value = 5.808e-06

# Because of that, lets look at a Kruskal test
kruskal.test(Cumulative_Respiration ~ site, data = resp) # Kruskal-Wallis chi-squared = 15.405, df = 2, p-value = 0.0004517

# Now let's visualize witha boxplot 
boxplot(Cumulative_Respiration ~ site, data = resp)

# Let's run a dunn test for respiration and site 
dunnTest(Cumulative_Respiration ~ site, data = resp) # still don't know if I should use a method or not 

# Need to do copy number pre and post thaw and maybe change in DOC or TDN? 

shapiro.test(cpnumb$copy_number_per_gsoil_PRE) #p-value = 0.005065, not normal
shapiro.test(resp4copy$Cumulative_Respiration) #p-value = 9.524e-06, not normal

kruskal.test(copy_number_per_gsoil_PRE ~ site, data = cpnumb) # chi-squared = 3.1581, df = 2, p-value = 0.2062

boxplot(copy_number_per_gsoil_PRE ~ site, data = cpnumb)

# Post thaw copy number analysis 
shapiro.test(cpnumb$copy_number_per_gsoil_POST) #p-value = 5.065e-06, not normal

kruskal.test(copy_number_per_gsoil_POST ~ site, data = cpnumb) # chi-squared = 6.7302, df = 2, p-value = 0.03456

boxplot(copy_number_per_gsoil_POST ~ site, data = cpnumb)


