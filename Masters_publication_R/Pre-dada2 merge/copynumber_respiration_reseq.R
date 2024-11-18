# Joy O'Brien
# Master's Research; Ernakovich Lab
# April 3, 2023


# Load the necessary libraries
library(vegan)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(e1071) #for skewness function
library(FSA)
library(devtools)
library(ggpubr)

# Read in the data from excel (change this to the repo when you move the data to repo)
resp4copy <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "Respiration_cum_forcopynumb")
cpnumb <- read_excel("~/GitHub/Masters_publication/incubation_physical_chemical.xlsx", sheet = "copy_numbers")

# Visualizing the data to start off 
scatter.smooth(cpnumb$copy_number_per_gsoil_PRE, y = resp4copy$Cumulative_Respiration)

# Checking for normality
shapiro.test(cpnumb$copy_number_per_gsoil_PRE) #p-value = 0.005065 Data is not normal
shapiro.test(resp4copy$Cumulative_Respiration) #p-value = 0.03419 data is not normal

# Since the data is not normally distributed we can run a KW test followed by a Dunn test
# Merge the data frames
respcpnumb <- merge(resp4copy,cpnumb) 

# For pre-thaw copy numbers and resp
kruskal.test(copy_number_per_gsoil_PRE ~ Cumulative_Respiration, data = respcpnumb) # Kruskal-Wallis chi-squared = 21, df = 21, p-value = 0.4599
# For pre-thaw copy numbers and site
kruskal.test(copy_number_per_gsoil_PRE ~ site, data = cpnumb) # not significant p = 0.2072

# For post thaw copy numbers and cumulative respiration
kruskal.test(copy_number_per_gsoil_POST ~ Cumulative_Respiration, data = respcpnumb) # Kruskal-Wallis chi-squared = 21, df = 21, p-value = 0.4589 SAME P VALUE?

# For post-thaw copy numbers and site
kruskal.test(copy_number_per_gsoil_POST ~ site, data = cpnumb) #significant p-value = 0.0156
dunn_test(copy_number_per_gsoil_POST ~ site, data = cpnumb)
kruskal.test(copy_number_per_gsoil_PRE ~ copy_number_per_gsoil_POST, data = cpnumb) # p-value = 0.4599

# Spearman Correlations
# Correlation for pre-thaw copy number and cumulative resp.
resp_pre_corr <- cor.test(x =respcpnumb$copy_number_per_gsoil_PRE, y = respcpnumb$Cumulative_Respiration, method = "spearman")
print(resp_pre_corr)

# Visualize here 
ggplot(respcpnumb, aes(x = copy_number_per_gsoil_PRE, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site)) +
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm")


#**************************************************************************
# FIGURE INCLUDING PRE AND POST THAW COPY NUMBER AND CUMULATIVE RESPIRATION
# Visualize both pre and post as a response to respiration (code provided by HHM): 
# Need to make the ns line dashed and the s line solid, need to flip the legend for the copy number type
# respcpnumb <- factor(respcpnumb, levels = rev(levels(respcpnumb)))
# MAKE SURE respcnumb is a df
pivot_longer(respcpnumb, cols = contains("copy_number"), names_to = "Pre_post_copynumber", values_to = "copy_number") %>% 
  ggplot(aes(respcpnumb, x = copy_number, y = Cumulative_Respiration, shape = Pre_post_copynumber, group = Pre_post_copynumber, color = site)) + 
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_x_log10() +
  theme_classic() +
  xlab("log10 Copy number g-1 soil") +
  ylab("Cumulative respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 12)) +
  scale_shape_discrete("Copy number type", labels = c("Post-thaw", "Pre-thaw")) + 
  scale_color_discrete("Site") +
  scale_linetype_manual(values = c("dashed"))

# CHat GPT
pivot_longer(respcpnumb, cols = contains("copy_number"), names_to = "Pre_post_copynumber", values_to = "copy_number") %>% 
  ggplot(aes(respcpnumb, x = copy_number, y = Cumulative_Respiration, shape = Pre_post_copynumber, group = Pre_post_copynumber, color = site)) + 
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  scale_x_log10() +
  theme_classic() +
  xlab("log10 Copy number g-1 soil") +
  ylab("Cumulative respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 12)) +
  scale_shape_discrete("Copy number type", labels = c("Post-thaw", "Pre-thaw")) + 
  scale_color_discrete("Site") +
  scale_linetype_manual(values = c("dashed"))
# Dashed insignifcant and solid significant
pivot_longer(respcpnumb, cols = contains("copy_number"), names_to = "Pre_post_copynumber", values_to = "copy_number") %>% 
  ggplot(aes(respcpnumb, x = copy_number, y = Cumulative_Respiration, shape = Pre_post_copynumber, group = Pre_post_copynumber, color = site, linetype = Pre_post_copynumber)) + 
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_log10() +
  theme_classic() +
  xlab("log10 Copy number g-1 soil") +
  ylab("Cumulative respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 12)) +
  scale_shape_discrete("Copy number type", labels = c("Post-thaw", "Pre-thaw")) + 
  scale_color_discrete("Site") +
  scale_linetype_manual(values = c("solid", "dashed"))

pivot_longer(respcpnumb, cols = contains("copy_number"), names_to = "Pre_post_copynumber", values_to = "copy_number") %>% 
  ggplot(aes(respcpnumb, x = copy_number, y = Cumulative_Respiration, shape = Pre_post_copynumber, group = Pre_post_copynumber, color = site, linetype = Pre_post_copynumber)) + 
  geom_point(size = 2.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_log10() +
  theme_classic() +
  xlab("log10 Copy number g-1 soil") +
  ylab("Cumulative respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 12)) +
  scale_shape_discrete("Copy number type", labels = c("Post-thaw", "Pre-thaw")) + 
  scale_color_discrete("Site") +
  scale_linetype_manual(values = c("solid", "dashed"), breaks = c("Pre-thaw", "Post-thaw"))

# most updated and recent code as of July 6, 2023
# Code that includes pre-thaw before post-thaw and the dashed/solid figure line
pivot_longer(respcpnumb, cols = contains("copy_number"), names_to = "Pre_post_copynumber", values_to = "copy_number") %>% 
  ggplot(aes(respcpnumb, x = copy_number, y = Cumulative_Respiration, shape = Pre_post_copynumber, group = Pre_post_copynumber, color = site, linetype = Pre_post_copynumber)) + 
  geom_point(size = 3.5) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_log10() +
  theme_classic() +
  xlab("log10 Copy number g-1 soil") +
  ylab("Cumulative respiration (µg C-CO2 g-1 dry soil)") +
  theme(text = element_text(size = 12)) +
  scale_shape_manual("Copy number type", labels = c("Post-thaw", "Pre-thaw"), values = c(2, 1)) +
  scale_color_discrete("Site") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  guides(shape = guide_legend(reverse = TRUE), linetype = "none") 

# Notes from Jessica regarding fixing this figure
# I like the colored sites and shapes for copy number. I don't love the colors for the lines. 
# If I recall, the pre-thaw is ns and the post thaw is significant. If so, let's make them both black lines, but make the n.s. line dashed.

# for the figure above need to move the order of pre-thaw and post-thaw and then make the 

http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
#*********************************************************************
# Correlation for post-thaw copy number and cumulative respiration
resp_post_corr <- cor.test(x = respcpnumb$copy_number_per_gsoil_POST, y = respcpnumb$Cumulative_Respiration, method = "spearman")
print(resp_post_corr)

# Visualize 
ggplot(respcpnumb, aes(x =copy_number_per_gsoil_POST, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm")

# FIGURE
# Cleaning up the figure for publication 
ggplot(respcpnumb, aes(x = copy_number_per_gsoil_POST, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site), size = 2.5) +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10() +
  xlab("log10 Post-thaw copy number g-1 soil") +
  ylab("log10 Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "black", se=TRUE) +
  theme(text = element_text(size = 12)) +
  scale_color_discrete("Site")

#
ggplot(respcpnumb, aes(x = copy_number_per_gsoil_PRE, y = Cumulative_Respiration)) + 
geom_point(aes(color = site), size = 2.5) +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10() +
  xlab("log10 Pre-thaw copy number g-1 soil") +
  ylab("log10 Cumulative Respiration (µg C-CO2 g-1 dry soil)") +
  geom_smooth(method = "lm", color = "black", se=TRUE, linetype = "dashed") +
  theme(text = element_text(size = 12)) +
  scale_color_discrete("Site")

# AND WE WILL RUN KRUSKALL ON BACTERIAL ABUNDANCE POST THAW BY SITE
kruskal.test(copy_number_per_gsoil_POST ~ site, data = respcpnumb) # significant: Kruskal-Wallis chi-squared = 6.7302, df = 2, p-value = 0.03456
dunnTest(copy_number_per_gsoil_POST ~ site, data = respcpnumb, method = 'bh') # not sure if I should be using bh or holms? 

# Now let's look to see if there is any significance pre and post thaw by site for copy numbers 
kruskal.test(copy_number_per_gsoil_PRE ~ copy_number_per_gsoil_POST, data = respcpnumb) # NOT SIGNIF Kruskal-Wallis chi-squared = 21, df = 21, p-value = 0.4589

# *************************************************************************************************************************************************************
# Creating the glm THIS NEEDS TO BE EDITED
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

# ANALYZING POST-THAW COPY NUMBERS AND RESPIRATION

# Visualizing the data to start off 
scatter.smooth(x = cpnumb$copy_number_per_gsoil_POST, y = resp4copy$Cumulative_Respiration)

# Checking for normality in the post thaw copy number values 
shapiro.test(cpnumb$copy_number_per_gsoil_POST) #p-value = 5.065e-06, not normal

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

# Run the glm, this needs to be edited
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

#*****************************************************************************

# Spearman correlation for copy numbers post thaw and copy numbers pre thaw

# First let's look at pre-thaw

# Make a QQ plot
qqnorm(x = cpnumb$copy_number_per_gsoil_PRE, y = resp4copy$Cumulative_Respiration, main = "Q-Q Plot")

# Run the corr test 
copyprecorr <- cor.test(x = cpnumb$copy_number_per_gsoil_PRE, y = resp4copy$Cumulative_Respiration, method = "spearman") #p value = 0.6012, rho = 0.1146245
print(copyprecorr)

precopyresp<- merge(cpnumb, resp4copy)

# Visualize
ggplot(precopyresp, aes(x = copy_number_per_gsoil_PRE, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, shape = core)) 
scale_x_log10() +
  scale_y_log10()

# Now we will look at copy number post thaw 

# Visualize a QQ plot first
qqnorm(x = copynumb$copy_number_per_gsoil_POST, y = resp4copy$Cumulative_Respiration, main = "Q-Q Plot")

copypostcorr <- cor.test(x = cpnumb$copy_number_per_gsoil_POST, y = resp4copy$Cumulative_Respiration, method = "spearman") # p-value = 0.001604, rho = 0.6304348
print(copypostcorr)

# Visualize
ggplot(precopyresp, aes(x = copy_number_per_gsoil_POST, y = Cumulative_Respiration)) + 
  geom_point(aes(color = site, shape = core)) +
  scale_x_log10() +
  scale_y_log10()

# LOOKING AT THE CUMULATIVE RESPIRATION BY ITSELF 
ggplot(data=precopyresp, aes(x=core, y=Cumulative_Respiration, group=1)) +
  geom_bar()

# END