library(tidyverse) # for data handling
library(purrr) # for pluck()
# To install propr:
# install.packages("devtools")
# devtools::install_github("tpq/propr")
library(propr) # for clr transformation

# function to remove 0s for transformation
transform_perc <- function(vec) {
  # Remove 0s for CLR transformation
  # See Smithson & Verkuilen 2006 (https://doi.org/10.1037/1082-989X.11.1.54)
  (vec * (length(vec) - 1) + 0.5) / length(vec)
}

ASV_table_fp <- "~/Downloads/ASV_table_prerare.csv"
ASV_table_fp <- "~/Downloads/ASV_table_prerare_specrownames.csv"
map_file_fp <- "~/Downloads/mapping_file.csv"


# read in data
ASV_table <- read_csv(file = ASV_table_fp)
mapping_table <- read_csv(file = map_file_fp)




# ASVS_of_interest
# mags_of_interest <- c("ASV_1", "ASV_7", "ASV_30") # can specify only a subset of the organisms if you don't need all of them for testing


abund_TM_RA_clr <- ASV_table %>% # a data frame with ASVs in rows and samples as columns; contains ASV ids as a column rather than rownames. The sample id column is called "SampleID"
  # translate dataframe to work with propr expectation
  column_to_rownames(var = "ASV_ID") %>%
  t() %>% as_tibble(rownames = "SampleID") %>% # flips the matrix so rows are samples and columns are ASVs. 
  column_to_rownames(var = "SampleID") %>% 
  # apply zeros transformation corrections - clr doesn't work if 0s are present
  mutate(across(everything(), transform_perc)) %>%
  # apply clr transformation
  propr(metric = "rho") %>%
  pluck("logratio") %>%
  rownames_to_column(var = "SampleID") %>%
  select(SampleID, everything())
  #select(SampleID, any_of(mags_of_interest)) # or can alternatively use: select(SampleID, everything()) if you don't have specific taxa you want to test

# put the transformed data into conventional formats:
# wide format
clr_transformed_wide_format <- abund_TM_RA_clr %>% column_to_rownames(var = "SampleID") %>%
  t() %>% as_tibble(rownames = "ASV_ID")
# long format
clr_transformed_long_format <- abund_TM_RA_clr %>%
  pivot_longer(-SampleID, names_to = "ASV_ID", values_to = "CLR_relabund")


# Write out data

write_csv(x = clr_transformed_wide_format, file = "~/Downloads/ASV_table_clr_transformed_wide_format.csv")
write_csv(x = clr_transformed_long_format, file = "~/Downloads/ASV_table_clr_transformed_long_format.csv")