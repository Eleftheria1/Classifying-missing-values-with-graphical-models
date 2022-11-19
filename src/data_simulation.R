######################### Simulate data ######################### 
# For each of the datasets to be simulated we follow the recipe below:
#   1 simulate the raw data
#   2 visualize
#   3 fit a regression
#   4 add missing values
#   5 visualize again

# Load the packages required:
library(tidyverse)
library(simstudy)

######################### Multivariate Normal ######################### 

# Simulate the raw data
simulate_multi_normal_dataset <- function(
    mu, sigma,
    target_formula,
    target_sd,
    n = 1000,
    seed = 2
) {
  var_names <- paste0("x", seq_along(mu))
  multi_normal_data_def <- defData(
    varname = "x1", dist = "normal", formula = mu[1], variance = sigma[1]
  )
  for (variable_count in 2:length(mu)) {
    multi_normal_data_def <- defData(
      multi_normal_data_def,
      varname = var_names[variable_count],
      dist = "normal",
      formula = mu[variable_count],
      variance = sigma[variable_count]
    )
  }
  multi_normal_data_def <- defData(
    multi_normal_data_def,
    varname = "y",
    formula = target_formula,
    dist = "normal",
    variance = target_sd
  )
  set.seed(seed)
  genData(n, multi_normal_data_def)
}

multi_normal <- simulate_multi_normal_dataset(
  mu = c(0, 5, 10, -1),
  sigma = c(1, 2, 3, 0.3),
  target_formula = "-3 * x1 + 5 * x2 - 0.5 * x3 + 5 * x4",
  target_sd = 3
)

# Visualize the raw data
GGally::ggpairs(
  multi_normal,
  columns = 2:ncol(multi_normal),
  lower = list(continuous = "density")
) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.ticks = element_blank(), 
    axis.text = element_blank()
  )

# Assess the fit of the raw data
lm(
  y ~ .,
  data = multi_normal %>% select(-id)
) %>%
  summary()


# Store the raw data in the final data list for the multivariate normal 
# simulations
multi_normal_data_list <- list(
  raw = list(
    data = list(multi_normal),
    missing = NA,
    missing_type = NA,
    relative_missingness = 0
  )
)

add_to_data_list <- function(
  data_list,
  simulation_name,
  data,
  missing,
  missing_type,
  relative_missingness
) {
  data_list[[simulation_name]] <- list(
    data = data,
    missing = missing,
    missing_type = missing_type,
    relative_missingness = relative_missingness
  )
  data_list
}
# Simulate missingness

relative_missingness <- c(0.1, 0.3, 0.6)


# MCAR for x2 for each relative missingness
simulation_list <- lapply(relative_missingness, function(rel_missingness) {
  missingness_x2 <- defMiss(varname = "x2", formula = rel_missingness)
  set.seed(2)
  miss_mat <- genMiss(
    multi_normal_data_list$raw$data[[1]],
    missingness_x2,
    idvars = "id"
  )
  genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")
})

# add simulated data to data list
multi_normal_data_list <- add_to_data_list(
  multi_normal_data_list,
  simulation_name = "mcar_x2",
  data = simulation_list,
  missing = "x2",
  missing_type = "mcar",
  relative_missingness = relative_missingness
)

save(
  multi_normal_data_list,
  file = paste0(here::here(), "/data/multi_normal.RData")
)





