# ANOVAs_site_script_reseq.R

# Joy O'Brien
# Master's work, Ernakovich lab
# March 28, 2023

# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(FSA) #installed for the dunntest
library(multcompView) # used for the tukey test
library(ggplot2) # to make pretty figs 
library(rstatix) 
library(multcompView) # allows us to obtain significance letters for post-hoc tests
library(rcompanion) # attempted to use this to get the significance codes from KW test
library(agricolae) # used to get the significance codes from the KW test

# Read in the data (this is all accessible in the incubation_physical_chemical.xlsx sheet on Ernakovich Lab Box Folder)
# Need to correct the file path below with the Git repo

cpnumb <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "copy_numbers")
resp <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum")
resp4copy <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Respiration_cum_forcopynumb")
resp_tctn <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "Resp_TC_TN")
ph_pre <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_pre")
ph_post <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_post")
ph_both <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "pH_both")
EC_pre <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_pre")
EC_post <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_post")
TC_TN <-  read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "TC_TN")
EC_both <-  read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "EC_both")
GWC_both <- read_excel("~/Desktop/incubation_physical_chemical.xlsx", sheet = "GWC_both")

# CUMULATIVE RESPIRATION

# Let's look at the differences in cumulative respiration across sites 
shapiro.test(resp$Cumulative_Respiration) #this is not normally distributed W = 0.89709, p-value = 0.01868

# Because of that, lets look at a Kruskal test
kruskal.test(Cumulative_Respiration ~ site, data = resp) # Kruskal-Wallis chi-squared = 16.805, df = 2, p-value = 0.0002243
kruskal(resp$Cumulative_Respiration, resp$site, group=TRUE, p.adj="bonferroni")$groups # use this to get significance codes for KW test
# CUMMULATIVE RESPIRATION IS STATISTICALLY SIGNIFICANT ACROSS SITES

# FIGURE
# Visualize cumulative respiration between sites with a box plot 
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
shapiro.test(resp4copy$Cumulative_Respiration) #p-value = 0.03419, not normal
# Kruskal on pre thaw copy number and site
kruskal.test(copy_number_per_gsoil_PRE ~ site, data = cpnumb) # chi-squared = 3.1581, df = 2, p-value = 0.2062

# Visualizing
boxplot(copy_number_per_gsoil_PRE ~ site, data = cpnumb)

# Post thaw copy number analysis 
shapiro.test(cpnumb$copy_number_per_gsoil_POST) #p-value = 5.065e-06, not normal
# Kruskal on post thaw copy number and site
kruskal.test(copy_number_per_gsoil_POST ~ site, data = cpnumb) # chi-squared = 8.321, df = 2, p-value = 0.0156
# Visualizing
boxplot(copy_number_per_gsoil_POST ~ site, data = cpnumb)

#********************************************************
#*ABIOTIC ANALYSES AND STATISTICS
#********************************************************
# pH pre-thaw ANOVAs

# Checking normality 
shapiro.test(ph_pre$pH_pre) # yay it's normally distributed, p-value = 0.3505
# Run ANOVA on pre-thaw pH
phpreaov <- aov(pH_pre ~ site, data = ph_pre)
print(phpreaov)
summary(phpreaov)

# Run TUKEY on pre-thaw pH anova 
phpretukey<- TukeyHSD(phpreaov, conf.level = 0.95)
print(phpretukey)
# Obtain significance codes for TUKEY
tukey_cld_pre_ph <-multcompLetters4(phpreaov, phpretukey)
print(tukey_cld_pre_ph)

# pH post-thaw ANOVAs

# Checking normality
shapiro.test(ph_post$pH_post) # p-value = 0.05193 WHEW
# Run ANOVA on post-thaw pH
phpostaov <- aov(pH_post ~ site, data = ph_post)
print(phpostaov)
summary(phpostaov)

# Run TUKEY
phposttukey<- TukeyHSD(phpostaov, conf.level = 0.95)
print(phposttukey)

# Comparing pre and post thaw pH
# ANOVA
phbothaov <- aov(pH_pre ~ pH_post, data = ph_both)
summary(phbothaov)
#***************************************************

# NOW WE WILL LOOK AT ELECTRICAL CONDUCTIVITY (EC) pre-thaw 

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

# Viewing results with significance codes
tukey_cld <-multcompLetters4(ecpostaov, ECposttukey)
print(tukey_cld)

#Comparing pre and post thaw EC via ANOVA 
# ANOVA
ecboth <- aov(EC_pre ~ EC_post, data = EC_both)
summary(ecboth)

#******************************************************
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
#******************************************************
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

#*****************************************************************
# Looking at C content and respiration 
shapiro.test(resp_tctn$Cumulative_Respiration) # p-value = 0.370
shapiro.test(resp_tctn$C) # p-value = 0.228
shapiro.test(resp_tctn$N) # p-value = 0.293

resp_c_aov <- aov(C ~ Cumulative_Respiration, data = resp_tctn)
summary(resp_c_aov)


# Looking at the relationship between C and Cumulative Resp 
resp_c_corr <- cor.test(x =resp_tctn$C, y = resp_tctn$Cumulative_Respiration, method = "spearman")
print(resp_c_corr)

# FIGURE: C and Cumulative Respiration regression
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
# by site
resp_n_aov <- aov(N ~ site, data = resp_tctn)
summary(resp_n_aov)
N_resptukey <- TukeyHSD(resp_n_aov)
print(N_resptukey)

# Looking at the relationship between N and Cumulative Resp 
resp_n_corr <- cor.test(x =resp_tctn$N, y = resp_tctn$Cumulative_Respiration, method = "spearman")
print(resp_n_corr)

# FIGURE: N and Cumulative Respiration Regression
# Making a pretty plot
ggplot(resp_tctn, aes(x = N, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site), size = 3.5) +
  theme_classic() +
  xlab("% N") +
  ylab("Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "dark gray", se = TRUE) +
  theme(text = element_text(size = 18)) +
  scale_color_discrete("Site") #changing the name of the legend title
#***********************************************************************
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
