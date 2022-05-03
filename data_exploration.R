# Joy O'Brien
# Visualizing physical/chemical Data 

#Install packages
#ggally for correlations, ggpairs

# Load packages 
library(readxl)
library (tidyverse)

# Loading data
data <-read_excel("~/Desktop/incubation_physical_chemical.xlsx")

#add a column to designate pre and post thaw samples
data2 <- data %>%
  mutate(pre_post_thaw = ifelse(test = grepl("PRE", Sample), no = "POST", yes = "PRE")) %>%
  mutate(site = ifelse(test=grepl("FL", Sample), yes= "FL",
                       no = ifelse(test=grepl("AT", Sample), yes="CRREL", no= "Utqiagvik")))
# see if there is a relationship between DOC and TDN (not statistically, just visually)
ggplot(data = data2, aes(x=DOC_ppm, y=TDN_ppm))+
  geom_point(aes(color=site, shape=pre_post_thaw))

ggplot(data=data2, aes(x=%_Moisture, y=pH))+
  geom_point(aes(color=site))
