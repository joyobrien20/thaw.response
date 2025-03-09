# Joy O'Brien
# Cumulative respiration figure
# October 4, 2023


library(ggplot2)
# Load data
data <- X20231018_CO2rates_incubation

# Subset data out by site 
CRREL_resp <- subset.data.frame(data, Site == "CRREL")

ggplot(data = data, aes(x = Day, y = Avg_cumresp, group = Site, color = Site)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  geom_errorbar(aes(ymin = Avg_cumresp - STD_dev_mean, ymax = Avg_cumresp + STD_dev_mean), width = 0.5, end = "both") +
  labs(y = "C-CO2, ug C-CO2/g dry soil", x = "Time, days") + 
  scale_color_manual(
  values = c("FL" = "green", "CRREL" = "red", "Ut" = "blue"), 
  labels = c("FL" = "Farmer's Loop", "Ut" = "Utqiagvik"))


  
#facet_wrap(~ Site, ncol = 1) 
# Change the line typ
# Basic line plot with points
ggplot(data=CO2_rates_FINAL_v2_copy$Site, aes(x = Day, y=C_CO2_rate, group=1)) +
  geom_line(aes(color = Site)) 
  #facet_wrap(~ Site, ncol = 1) 
# Change the line type
ggplot(data=df, aes(x=dose, y=len, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()
# Change the color
ggplot(data=df, aes(x=dose, y=len, group=1)) +
  geom_line(color="red")+
  geom_point()
