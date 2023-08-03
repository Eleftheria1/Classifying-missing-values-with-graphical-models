###############################################################################
# This r- file uses the R-code from imputation_utils.R to perform the imputations
###############################################################################
#load required packages
library(Amelia)
library(clarify)
library(tidyverse)
library(mice)
library(bestNormalize)

#load required data and src code
base_datasets <- c("mixed_data.RData", "multi_normal.RData")
all_data_files <- list.files(path = "data")
experiment_results <- all_data_files[!(all_data_files %in% base_datasets)]

load(paste0("data/", base_datasets[1]))
load(paste0("data/", base_datasets[2]))
source("src/imputation_utils.R")

set.seed(123)
 #imputed_mixed_data_list <- overall_impute(
 #  data_list = mixed_data_list,
 #  nominals = c("x6", "x7", "x8", "x9"),
 #  sqrts = c("x3"),
 #  m = 500
 #) 
#save(imputed_mixed_data_list, file = "data/imputed_mixed_data.RData")
load("data/imputed_mixed_data.RData")

# imputed_normal_data_list <- overall_impute(
#   data_list = multi_normal_data_list
# )
#save(imputed_normal_data_list, file = "data/imputed_normal_data.RData")
load("data/imputed_normal_data.RData")


# Optional Normalization -------------------------------------------------
trans_cols <- c(
  "x2",
  "x6",
  "x7",
  "x8",
  "x9"
)
set.seed(222)
normalized_mixed_data_list <- normalize_columns(
  mixed_data_list, columns = trans_cols
)

imputed_norm_mixed_data_list <- overall_impute(
  data_list = normalized_mixed_data_list,
  nominals = c("x6", "x7", "x8", "x9"),
  sqrts = c("x3"),
  m = 500
)


imputed_norm_mixed_data_list <- retransform_columns(
  imputed_norm_mixed_data_list, columns = trans_cols
)
save(imputed_norm_mixed_data_list, file = "data/imputed_norm_mixed_data.RData")
