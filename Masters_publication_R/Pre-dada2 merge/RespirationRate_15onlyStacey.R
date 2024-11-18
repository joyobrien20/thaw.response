# This script will import a csv file of CO2 respiration rates and calculate 
# average respiration rates, Q10, and ANOVA (need to add this)

####### all data #######
# First change the working directory to where the data is.
setwd("/Users/staceydoherty/Box/Lab logistics/Doherty Thesis/Thesis/Incubation study/R files")

# Import data from .csv file. Data file should have one header row and columns containing each categorical and response variable.
# Column 27 was deleted for the 4C data because that was the day the incubator was turned off.
rates = read.csv("CO2_rates_forR_15C_only.csv", header=TRUE, sep = ",", dec = ".", na.strings = "NA")

# This may be unnecessary, but the above data was converted to a matrix.
mat_rates = as.data.frame(rates)

# A new matrix was created with the mean and sd for each sample. 
mat_rates$avg_rate = rowMeans(mat_rates[,c(6:37)], na.rm = TRUE, dims = 1)

# This appended the sd to the existing matrix.
mat_rates$sd = apply(mat_rates[,c(6:37)], 1, sd, na.rm = TRUE)

# Standard error was calculated and appended to the matrix.
mat_rates$se = mat_rates$sd/sqrt(ncol(mat_rates[,c(6:37)])-rowSums(is.na(mat_rates)))

# Plot the mean and se for each jar to see how rates fluctuated over time
# Load ggplot2
library(ggplot2)
# Values were plotted with error bars. 
ggplot(mat_rates, aes(x=sample_ID, y=avg_rate, group=1)) + geom_point() + geom_errorbar(data=mat_rates, 
  aes(ymin=avg_rate-se, ymax=avg_rate+se), width=0.2) + theme_classic()
# Error bars look small except when rates were large, I think its okay to go ahead and average by depth!

# Rates were plotted as a box plot averaged by depth and temp
ggplot(mat_rates, aes(x=depth, y=avg_rate, fill=temp)) + geom_boxplot() + theme_classic() + scale_fill_manual("temp", values = c("15C" = "white", "4C" = "grey"))

######### Also created a box plot where only days 22-193 were included ######
# since that is where rates visually begin to level off.
# To do this run the above code but change the range for calcs to c(,15:37)
# These are the data used in Q10 calcs

# make new data frame with days 22-193
# remove days 1-21.
mat_rates_stable = as.data.frame(rates[, -c(6:14)])

# A new matrix was created with the mean and sd for each sample. 
mat_rates_stable$avg_rate = rowMeans(mat_rates[,c(6:27)], na.rm = TRUE, dims = 1)

# This appended the sd to the existing matrix.
mat_rates_stable$sd = apply(mat_rates[,c(6:27)], 1, sd, na.rm = TRUE)

# Standard error was calculated and appended to the matrix.
mat_rates_stable$se = mat_rates$sd/sqrt(ncol(mat_rates[,c(6:27)])-rowSums(is.na(mat_rates)))

# Plot the mean and se for each jar to see how rates fluctuated over time
# Load ggplot2
library(ggplot2)
# Values were plotted with error bars. 
ggplot(mat_rates_stable, aes(x=sample_ID, y=avg_rate, group=1)) + geom_point() + geom_errorbar(data=mat_rates_stable, aes(ymin=avg_rate-se, ymax=avg_rate+se), width=0.2) + theme_classic()
# Error bars look small except when rates were large, I think its okay to go ahead and average by depth!

# Rates were plotted as a box plot averaged by depth and temp
ggplot(mat_rates_stable, aes(x=depth, y=avg_rate, fill=temp)) + geom_boxplot() + theme_classic() + scale_fill_manual("temp", values = c("15C" = "white", "4C" = "grey"))

#### add 2 way anova here!

######## Compute the analysis of variance - one way by depth
res2.aov <- aov(avg_rate ~ depth, data = mat_rates_stable)

## Summary of the analysis
summary(res2.aov)

######### Check Assumptions of the test

## Assumption 1: Normality of the residual
plot(res2.aov,2) #Q-Q plot

shapiro.test(x = residuals(object = res2.aov)) # Shapiro-Wilk test

## Assumption 2: Homogeineity of Variance
library(car) # Use Levene test from the car library

leveneTest(avg_rate ~ temp*depth, data = mat_rates_stable)

# since the normality assumption failed, we need to do Kruskal wallis test

####################### Non-parametric laternative to two-ways ANOVA
## This isn to be used when at least one of the assumptions is not satisfied.

kruskal.test(avg_rate ~ temp, data = mat_rates_stable) # Kruskal-Wallis rank test
kruskal.test(avg_rate ~ depth, data = mat_rates_stable) # Kruskal-Wallis rank test


##### Q10 cals #####

# coudn't figure out calculation so I did them in excel and imported the Q10 data table to do stats

# First change the working directory to where the data is.
setwd("/Users/staceydoherty/Box/Lab logistics/Doherty Thesis/Thesis/Incubation study/R files")

# Import data from .csv file. Data file should have one header row and columns containing each categorical and response variable.
q10 = read.csv("q10_forR.csv", header=TRUE, sep = ",", dec = ".")

# This may be unnecessary, but the above data was converted to a matrix.
q10 = as.data.frame(q10)

# q10 values were plotted as a box plot averaged by depth
ggplot(q10, aes(x=depth, y=Q10)) + geom_boxplot() + theme_classic()


######## Compute the analysis of variance - this is the one way anova by depth!
q10.aov <- aov(Q10 ~ depth, data = q10)

## Summary of the analysis
summary(q10.aov)


######### Check Assumptions of the test

## Assumption 1: Normality of the residual
plot(q10.aov,2) #Q-Q plot

shapiro.test(x = residuals(object = q10.aov)) # Shapiro-Wilk test

## Assumption 2: Homogeineity of Variance
library(car) # Use Levene test from the car library

leveneTest(Q10 ~ depth, data = q10)

####################### Tukey Multiple comparison test

## Note: Multiple comparison is only used when you have significance (i.e. P<0.05)

TUKEY_q10 <- TukeyHSD(x = q10.aov, 'depth', conf.level=0.95)

TUKEY_q10

## A better visualization of multiple comparison results using letters

library(multcomp) ## library needed for the function
library(multcompView) ## library needed for the function

## Create a function to generate letters according to Tukey results
generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

# visualize it
generate_label_df(TUKEY_q10,'depth')



