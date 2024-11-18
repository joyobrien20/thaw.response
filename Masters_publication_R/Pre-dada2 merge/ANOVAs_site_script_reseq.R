# ANOVAs_site_script_reseq.R
#******EDIT!!!!*******
# Joy O'Brien
# Master's work, Ernakovich lab
# May 4, 2022

# THIS IS THE ANOVAS FOR SITE COMPARISON SCIPT
# 
# Data we are working with here: We will be running ANOVAs on DOC data post thaw, TDN data post thaw, 
# copy numbers pre thaw and copy numbers post thaw and cum respiration across sites

# The reason we are NOT including pre thaw data here is because the pre thaw data is not properly replicated (due to lack of soil to replicate). 

#* Random script notes*
#* 
# Visualizing/checking for theoretical model distribution (https://www.statology.org/anova-assumptions/)
# hist(doctdn$DOC_post)
# qqnorm(anova_docpost$residuals)


# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(FSA) #installed for the dunntest
library(multcompView) # used for the tukey test
library(ggplot2) # to make pretty figs 
library(rstatix) 
library(multcompView)
library(rcompanion) # attempted to use this to get the significance codes from KW test
library(agricolae) # used to get the significance codes from the KW test

# Read in the data (this is all accessible in the incubation_physical_chemical.xlsx sheet on Ernakovich Lab Box Folder)
doctdn <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "DOC_TDN")
resp <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "Respiration_cum")
ph_pre <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "pH_pre")
ph_post <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "pH_post")
ph_both <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "pH_both")
EC_pre <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "EC_pre")
EC_post <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "EC_post")
TC_TN <-  read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "TC_TN")
EC_both <-  read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "EC_both")
GWC_both <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "GWC_both")
texture <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "Texture")
doctdn_diff <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "DOC_TDN_diff")
resp_tctn <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "Resp_TC_TN")

# FIRST: RUNNING STATS ON THE DIFFERENCE IN DOC AND TDN PRE AND POST THAW

# DOC
# Checking for normality with a Shapiro test for DOC diff
shapiro.test(doctdn_diff$DOC_diff) #p-value = 0.0056; this data is not NORMAL

# The data is NOT normally distributed so I am going to try a Kruskall-Wallis test 
kruskal.test(DOC_diff ~ site, data = doctdn_diff) # Kruskal-Wallis chi-squared =  chi-squared = 4.16, df = 2, p-value = 0.1249
# The difference in DOC pre and post thaw is not significant across sites 

# May want to do a pairwise wilcox test for this ? http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
dunnTest(DOC_post ~ site, data = doctdn_diff, method = "bonferroni") # Check to see if the method is necessary? If it needs to be changed 

#let's visualize DOC_diff content across site 
boxplot(DOC_diff ~ site, data = doctdn_diff)

# For thesis and pub: making the doc_diff boxplot prettier 

doc <- ggplot(doctdn_diff, aes(x = site, y =DOC_diff)) + 
  geom_boxplot() +
  theme_classic() +
  xlab("Site") +
  ylab("Change in DOC (ppm)") +
  theme(text = element_text(size = 18))

plot(doc)
# TDN
# Checking for normality with a Shapiro test for TDN diff 
shapiro.test(doctdn_diff$TDN_diff) # p-value = 0.01293 # this is not normal as expected

# Therefore, we will run a Kruskall-Wallis test 
kruskal.test(TDN_diff ~ site, data = doctdn_diff) # Kruskal-Wallis chi-squared = 10.445, df = 2, p-value = 0.005394

kruskal(doctdn_diff$TDN_diff, doctdn_diff$site, group=TRUE, p.adj="bonferroni")$groups # use this to get the letter significance codes
# The difference in TDN pre and post thaw is statistically significant across sites

# Since the p-value of the Kruskal for TDN_diff is significant, we need to run a dunntest. 

# Now let's look at a dunn test
tdndunn <- dunnTest(TDN_diff ~ site, data = doctdn_diff, method = "bh")
print(tdndunn)
k_test <- k$doctdn_diff
cldList(comparison = k_test$Comparison,
        p.value    = PT$P.adj,
        threshold  = 0.05)

# Now let's visualize with a box plot
boxplot(TDN_diff ~ site, data = doctdn_diff)

# Making it better for pub and thesis 
tdn <- ggplot(doctdn_diff, aes(x = site, y = TDN_diff)) + 
  geom_boxplot() +
  theme_classic() +
  xlab("Site") +
  ylab("Change in TDN (ppm)") +
  theme(text = element_text(size = 18))
plot(tdn)
# DOC difference box plot 
geom_boxplot(doctdn_diff)

# Comparing TDN and DOC together 
kruskal.test(DOC_diff ~ TDN_diff, data = doctdn_diff)
#**********************************************************************

# NOW TRANSITIONING TO CUMULATIVE RESPIRATION

# Let's look at the differences in cumulative respiration across sites 
shapiro.test(resp$Cumulative_Respiration) #this is not normally distributed W = 0.89709, p-value = 0.01868

# Because of that, lets look at a Kruskal test
kruskal.test(Cumulative_Respiration ~ site, data = resp) # Kruskal-Wallis chi-squared = 16.805, df = 2, p-value = 0.0002243
kruskal(resp$Cumulative_Respiration, resp$site, group=TRUE, p.adj="bonferroni")$groups # use this to get significance codes for KW test
# CUMMULATIVE RESPIRATION IS STATISTICALLY SIGNIFICANT ACROSS SITES

# Now let's visualize with a box plot 
boxplot(Cumulative_Respiration ~ site, data = resp)

ggplot(resp, aes(x = site, y = Cumulative_Respiration)) + 
  geom_boxplot() +
  theme_classic() +
  xlab("Site") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 18))

# Let's run a dunn test for respiration and site 
dunnTest(Cumulative_Respiration ~ site, data = resp, method ='bh') 

# Checking to see if the copy number data pre thaw is normally distributed
shapiro.test(cpnumb$copy_number_per_gsoil_PRE) #p-value = 0.005065, not normal
# Checking to see if the cumulative respiration data is normally distributed
shapiro.test(resp4copy$Cumulative_Respiration) #p-value = 9.524e-06, not normal
# Kruskal on pre thaw copy number and site
kruskal.test(copy_number_per_gsoil_PRE ~ site, data = cpnumb) # chi-squared = 3.1581, df = 2, p-value = 0.2062
# Visualizing
boxplot(copy_number_per_gsoil_PRE ~ site, data = cpnumb)

# Post thaw copy number analysis 
shapiro.test(cpnumb$copy_number_per_gsoil_POST) #p-value = 5.065e-06, not normal
# Kruskal on post thaw copy number and site
kruskal.test(copy_number_per_gsoil_POST ~ site, data = cpnumb) # chi-squared = 6.7302, df = 2, p-value = 0.03456
# Visualizing
boxplot(copy_number_per_gsoil_POST ~ site, data = cpnumb)
#********************************************************
#*ABIOTIC ANALYSES AND STATISTICS
#********************************************************
# pH pre-thaw ANOVA

# Checking normality 
shapiro.test(ph_pre$pH_pre) # yay it's normally distributed, p-value = 0.3505
# Run ANOVA
phpreaov <- aov(pH_pre ~ site, data = ph_pre)
print(phpreaov)
summary(phpreaov)
# Run TUKEY
phpretukey<- TukeyHSD(phpreaov, conf.level = 0.95)
print(phpretukey)
# Obtain significance codes for TUKEY
tukey_cld_pre_ph <-multcompLetters4(phpreaov, phpretukey)
print(tukey_cld_pre_ph)
# $site
#diff        lwr          upr     p adj
#FL-CRREL        -0.155 -0.5992086  0.289208596 0.4216976
#Utqiagvik-CRREL -0.590 -1.0342086 -0.145791404 0.0233858
#Utqiagvik-FL    -0.435 -0.8792086  0.009208596 0.0527796

# pH post-thaw ANOVA 

# Checking normality
shapiro.test(ph_post$pH_post) # p-value = 0.05193 WHEW
# ANOVA
phpostaov <- aov(pH_post ~ site, data = ph_post)
print(phpostaov)
summary(phpostaov)

# Run TUKEY
phposttukey<- TukeyHSD(phpostaov, conf.level = 0.95)
print(phposttukey)
#$site
#diff        lwr         upr     p adj
#FL-CRREL        -0.140 -0.3624316  0.08243164 0.1510352
#Utqiagvik-CRREL -0.855 -1.0774316 -0.63256836 0.0011239
#Utqiagvik-FL    -0.715 -0.9374316 -0.49256836 0.0018593

# Comparing pre and post thaw pH
# ANOVA
phbothaov <- aov(pH_pre ~ pH_post, data = ph_both)
summary(phbothaov)
# TUKEY
phbothtukey <- TukeyHSD(phbothaov, conf.level = 0.95) #DO I ACTUALLY HAVE TO RUN THIS? COME BACK TO THIS
print(phbothtukey)

# NOW WE WILL LOOK AT ELECTRICAL CONDUCTIVITY (EC)
# EC pre-thaw ANOVA

# Checking for normality
shapiro.test(EC_pre$EC_pre) # p-value = 0.8672 
# ANOVA
ecpreaov <- aov(EC_pre ~ site, data = EC_pre)
summary(ecpreaov)
# Not significant so we don't need to run a Tukey

# EC post-thaw ANOVA 

# Checking for normality
shapiro.test(EC_post$EC_post) # p-value = 0.4292
# ANOVA
ecpostaov <- aov(EC_post ~ site, data = EC_post)
summary(ecpostaov)
# TUKEY
ECposttukey <- TukeyHSD(ecpostaov, conf.level = 0.95)
print(ECposttukey)

# this library allows us to view the significance codes associated with TUKEY
library(multcompView)
# Viewing results with significance codes
tukey_cld <-multcompLetters4(ecpostaov, ECposttukey)
print(tukey_cld)

#Comparing pre and post thaw EC via ANOVA 
# ANOVA
ecboth <- aov(EC_pre ~ EC_post, data = EC_both)
summary (ecboth)
# TUKEY
ECbothtukey <- TukeyHSD(ecboth, conf.level = 0.95)
print(ECposttukey)

# Comparing % C across site
# Checking for normality
shapiro.test(TC_TN$C) # p = 0.2355
TCaov <- aov(C ~ site, data = TC_TN)
summary(TCaov)
# Tukey
TCtukey <- TukeyHSD(TCaov, conf.level = 0.95)
print(TCtukey)

# obtaining results with significance codes
tukey_tc_tld <-multcompLetters4(TCaov, TCtukey)
print(tukey_tc_tld)

# Comparing % N across site 
# Checking for normality
shapiro.test(TC_TN$N) # p = 0. 3077 
# ANOVA
TNaov <- aov(N ~ site, data = TC_TN)
summary(TNaov)
# Tukey
TNtukey <- TukeyHSD(TNaov, conf.level = 0.95)
print(TNtukey)
# obtaining results with significance codes
tukey_tn <-multcompLetters4(TNaov, TNtukey)
print(tukey_tn)

# Now we are going to look at GWC across sites pre and post thaw 

# Pre-thaw GWC
# Checking for normality
shapiro.test(GWC_both$GWC_pre) # p-value = 0.6917
# ANOVA
GWC_preaov <- aov(GWC_pre ~ site, data = GWC_both)
summary(GWC_preaov) # Not significant across sites, so I will not run a Tukey

# Post-thaw GWC
# normality
shapiro.test(GWC_both$GWC_post) # p-value = 0.4414
# ANOVA
GWC_postaov <- aov(GWC_post~ site, data = GWC_both)
summary(GWC_postaov) 

# Comparing pre and post thaw
# anova
GWCchange <- aov(GWC_pre ~ GWC_post, data = GWC_both)
summary(GWCchange)
# tukey
gwctukey <- tukey_hsd(GWCchange)
print(GWCtukey)


# Anova on percent sand, silt, clay * NOTE: THIS DATA IS NOT USED IN THESIS/PUBLICATION
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

# Looking at C content and respiration 
shapiro.test(resp_tctn$Cumulative_Respiration) # p-value = 0.370
shapiro.test(resp_tctn$C) # p-value = 0.228
shapiro.test(resp_tctn$N) # p-value = 0.293

resp_c_aov <- aov(C ~ Cumulative_Respiration, data = resp_tctn)
summary(resp_c_aov)
Cresptukey <- TukeyHSD(resp_c_aov)
print(Cresptukey)

# Looking at the relationship between C and Cumulative Resp 
resp_c_corr <- cor.test(x =resp_tctn$C, y = resp_tctn$Cumulative_Respiration, method = "spearman")
print(resp_c_corr)

# Making a pretty plot
ggplot(resp_tctn, aes(x = C, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site)) +
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm")

ggplot(resp_tctn, aes(x = C, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site), size = 3.5) +
  theme_classic() +
  xlab("% C") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 18)) +
  scale_color_discrete("Site") #changing the name of the legend title

# Looking at the relationship between N and cumulative respiration 
resp_n_aov <- aov(N ~ site, data = resp_tctn)
summary(resp_n_aov)
N_resptukey <- TukeyHSD(resp_n_aov)
print(N_resptukey)

# Looking at the relationship between N and Cumulative Resp 
resp_n_corr <- cor.test(x =resp_tctn$N, y = resp_tctn$Cumulative_Respiration, method = "spearman")
print(resp_n_corr)

# Making a pretty plot
ggplot(resp_tctn, aes(x = N, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site), size = 3.5) +
  theme_classic() +
  xlab("% N") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 18)) +
  scale_color_discrete("Site") #changing the name of the legend title

# Trying GLM since this data is normally distributed 
# C
linearmodel_Cresp <- lm(Cumulative_Respiration ~ C, data = resp_tctn)

print(linearmodel_Cresp)
summary(linearmodel_Cresp)
# cor.test(resp_tctn$C, resp_tctn$Cumulative_Respiration, method = "pearson")

ggplot(linearmodel_Cresp, aes(x = C, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, size = 3.5)) +
  theme_classic() +
  xlab("% C") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 18)) +
  scale_color_discrete("Site") #changing the name of the legend title

# N
linearmodel_Nresp <- lm(Cumulative_Respiration ~ N, data = resp_tctn)
summary(linearmodel_Nresp)
print(linearmodel_Nresp)
lm.beta(linearmodel_Nresp)


shapiro.test(resp_tctn$N)
Nanova <- aov(N ~ Cumulative_Respiration, data = resp_tctn)
summary(Nanova)

ggplot(linearmodel_Nresp, aes(x = N, y = Cumulative_Respiration)) + 
  geom_point(aes(color = linearmodel_Cresp[["data"]][["site"]]), size = 2.5) +
  theme_classic() +
  xlab("% N") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 12)) +
  scale_color_discrete("Site") #changing the name of the legend title


# Checking the relationship using all of the replicates (Jessica asked to make sure that this relationship is strong)
respcheck <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Resp_TC_TN_check")

shapiro.test(respcheck$Cumulative_Respiration) # not normally distributed
resp_C_corr_test <- cor.test(x =respcheck$C, y = respcheck$Cumulative_Respiration, method = "spearman")
print(resp_C_corr_test)


ggplot(respcheck, aes(x = C, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site), size = 2.5) +
  theme_classic() +
  xlab("% C") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 12)) +
  scale_color_discrete("Site") #changing the name of the legend title

linearmodel_Crespcheck <- lm(Cumulative_Respiration ~ C, data = respcheck)
print(linearmodel_Crespcheck)
summary(linearmodel_Crespcheck)

# END
