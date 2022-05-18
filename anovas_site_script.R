# Joy O'Brien
# Master's work
# May 4, 2022

# THIS IS THE ANOVAS FOR SITE COMPARISON SCIPT
# 
# Data we are working with here: We will be running ANOVAs on DOC data post thaw, TDN data post thaw, copy numbers pre thaw and copy numbers post thaw and cum respiration across sites

# The reason we are NOT including pre thaw data here is because the pre thaw data is not properly replicated (due to lack of soil to replicate). 
# RUN A TUKEY TEST AFTER THE ANOVA TEST, NON PARAMETRIC ALTERNATIVE: SHAPIRO, KRUSKALL, DUNN 

# TO DO:
# NEED TO CHECK TO SEE IF YOU ARE RUNNING THE KRUSKALL WITH THE RIGHT PARAMETERS

# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(FSA) #installed for the dunntest
library(multcompView) # used for the tukey test

# Read in the data 
doctdn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "DOC_TDN")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum (2)")
ph_pre <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_pre")
ph_post <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_post")
ph_both <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_both")
EC_pre <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_pre")
EC_post <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_post")
TC_TN <-  read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "TC_TN")
EC_both <-  read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_both")
GWC_both <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "GWC_both")
texture <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Texture")
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
#**********************************************************
#*ABIOTIC ANALYSES AND STATISTICS
#**********************************************************
# pH pre-thaw ANOVA
shapiro.test(ph_pre$pH_pre) # yay it's normally distributed, p-value = 0.3505
phpreaov <- aov(pH_pre ~ site, data = ph_pre)
print(phpreaov)
summary(phpreaov)
phpretukey<- TukeyHSD(phpreaov, conf.level = 0.95)
print(phpretukey)

# $site
#diff        lwr          upr     p adj
#FL-CRREL        -0.155 -0.5992086  0.289208596 0.4216976
#Utqiagvik-CRREL -0.590 -1.0342086 -0.145791404 0.0233858
#Utqiagvik-FL    -0.435 -0.8792086  0.009208596 0.0527796



# pH post-thaw ANOVA 
shapiro.test(ph_post$pH_post) # p-value = 0.05193 WHEW
phpostaov <- aov(pH_post ~ site, data = ph_post)
print(phpostaov)
plot(phpretukey, las=1)
summary(phpostaov)
phposttukey<- TukeyHSD(phpostaov, conf.level = 0.95)
print(phposttukey)
  #$site
  #diff        lwr         upr     p adj
  #FL-CRREL        -0.140 -0.3624316  0.08243164 0.1510352
  #Utqiagvik-CRREL -0.855 -1.0774316 -0.63256836 0.0011239
  #Utqiagvik-FL    -0.715 -0.9374316 -0.49256836 0.0018593
  
# Comparing pre and post thaw pH
phbothaov <- aov(pH_pre ~ pH_post, data = ph_both)
summary(phbothaov)
phbothtukey <- TukeyHSD(phbothaov, conf.level = 0.95) #DO I ACTUALLY HAVE TO RUN THIS? COME BACK TO THIS
print(phbothtukey)

# EC pre-thaw ANOVA
shapiro.test(EC_pre$EC_pre) # p-value = 0.8672
ecpreaov <- aov(EC_pre ~ site, data = EC_pre)
summary(ecpreaov)
# Not significant so we don't need to run a Tukey

# EC post-thaw ANOVA 
shapiro.test(EC_post$EC_post) # p-value = 0.4292
ecpostaov <- aov(EC_post ~ site, data = EC_post)
summary(ecpostaov)
ECposttukey <- TukeyHSD(ecpostaov, conf.level = 0.95)
print(ECposttukey)

#Comparing pre and post thaw EC via ANOVA 
ecboth <- aov(EC_pre ~ EC_post, data = EC_both)
summary (ecboth)
ECbothtukey <- TukeyHSD(ecboth, conf.level = 0.95)
print(ECposttukey)

# Comparing % C across site 
shapiro.test(TC_TN$`%C`) # p = 0.2355
TCaov <- aov(`%C` ~ site, data = TC_TN)
summary(TCaov)
TCtukey <- TukeyHSD(TCaov, conf.level = 0.95)
print(TCtukey)

# Comparing % N across site 
shapiro.test(TC_TN$`%N`) # p = 0. 3077 
TNaov <- aov(`%N` ~ site, data = TC_TN)
summary(TNaov)
TNtukey <- TukeyHSD(TNaov, conf.level = 0.95)
print(TNtukey)

# Now we are going to look at GWC across sites pre and post thaw 

# Pre-thaw GWC
shapiro.test(GWC_both$GWC_pre) # p-value = 0.6917
GWC_preaov <- aov(GWC_pre ~ site, data = GWC_both)
summary(GWC_preaov) # Not significant across sites, so I will not run a Tukey

# Post-thaw GWC
shapiro.test(GWC_both$GWC_post) # p-value = 0.4414
GWC_postaov <- aov(GWC_post~ site, data = GWC_both)
summary(GWC_postaov) 

# Anova on percent sand, silt, clay 
# Clay
shapiro.test(texture$Clay) #p-value = 0.2864
clayaov <- aov(Clay ~ site, data = texture)
summary(clayaov) #p-value = 0.128 not statistically signif. 

# Sand
shapiro.test(texture$Sand) #p-value = 0.9378
sandaov <- aov(Sand ~ site, data = texture)
summary(sandaov) #p-value = 0.0594, not statistically significant 

# Silt 
shapiro.test(texture$Silt) #p-value = 0.1183
siltaov <- aov(Silt ~ site, data = texture)
summary(siltaov) # p-value = 0185, not statistically signif 



