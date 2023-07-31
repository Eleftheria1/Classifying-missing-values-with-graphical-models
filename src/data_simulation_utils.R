# Load the packages required:
library(tidyverse)
library(ggExtra)
library(latex2exp)
library(simstudy)

# Global Params
relative_missingness <- c(0.1, 0.3, 0.6)
# Normal Dataset
mu = c(0, 5, 10, -1, 3)
sigma = c(1, 2, 3, 0.3, 1)
target_formula = "-3 * x1 + 5 * x2 - x3 + 5 * x4 + 2 * x5"
target_sd = 3

# Generic Utils

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

#' Add missingness indicators to datasets of the data list
#'
#' @param data_list 
#'
#' @return amended data list
create_miss_ind <- function(data_list) {
  for (i in 1:length(data_list)) {
    for (k in 1:length(data_list[[i]]$data)) {
      df <- data_list[[i]]$data[[k]]
      missing_column_indicator <- colSums(is.na(df)) > 0
      missing_columns <- names(missing_column_indicator)[missing_column_indicator]
      for (missing_column in missing_columns) {
        df[[paste0("missing_", missing_column)]] <- if_else(
          is.na(df[[missing_column]]), 1, 0
        )
      }
      data_list[[i]]$data[[k]] <- df
    }
  }
  data_list
}

create_missingness_indicator <- function(
    raw_df, from, rel_missingness,
    seed = 2, sig = TRUE, multiplier = 1.5
) {
  set.seed(seed)
  if (sig) {
    sig <- function(x) { 1 / (1 + exp(-x)) }
    studentized_from <- (raw_df[[from]] - mean(raw_df[[from]])) / sd(raw_df[[from]])
    sample_probs <- as.numeric(
      sig(multiplier * studentized_from)
    )
  } else {
    ecdf_from <- ecdf(raw_df[[from]])
    sample_probs <- as.numeric(
      ecdf_from(raw_df[[from]])
    )
  }
  
  missing_indicator <- sample(
    seq_along(sample_probs),
    size = round(rel_missingness * length(sample_probs)),
    replace = FALSE,
    prob = sample_probs
  )
  seq_along(sample_probs) %in% missing_indicator
}

#Change names of variables to Latex code
var2TeX <- function(string, tex = TRUE) {
  if (string == "y") {
    return(string)
  }
  raw <- paste0(
    "$",
    str_extract(string, "([a-z]*)"),
    "_",
    str_extract(string, "([0-9]+)"),
    "$"
  )
  if (!tex) {
    return(raw)
  }
  TeX(raw)
}

# plot_single_missing: Generates a plot to explore missingness patterns in a dataset.
# Parameters:
#   - data_list: A list containing simulation data
#   - simulation_name: The name of the simulation within the data_list to be plotted.
#   - relative_missingness: The relative missingness level to be visualized.
#   - comparison_variable: The variable for comparison on the y-axis of the plot.
#   - alpha_non_missing: Alpha value (transparency) for non-missing data points (default = 1).
#   - density: Boolean; if TRUE, a 2D density plot is created, otherwise, points are plotted.
#   - missingness_index: Index to select a specific missingness variable (default = 1).
#   - marginal: Boolean; if TRUE and density is FALSE, marginal density plots are added.
#   - fixed_coords: Boolean; if TRUE, forces the plot's aspect ratio to be 1:1.

plot_single_missing <- function(
    data_list,
    simulation_name, 
    relative_missingness,
    comparison_variable = "x1",
    alpha_non_missing = 1,
    density = FALSE,
    missingness_index = 1,
    marginal = FALSE,
    fixed_coords = FALSE
) {
  rel_missingness_index <- which(
    data_list[[simulation_name]]$relative_missingness == relative_missingness
  )
  missing_variable <- data_list[[simulation_name]]$missing
  if (length(missing_variable) > 1) {
    missing_variable <- missing_variable[missingness_index]
  }
  mean_x <- mean(data_list$raw$data[[1]][[missing_variable]])
  mean_y <- mean(data_list$raw$data[[1]][[comparison_variable]])
  plot <- data_list$raw$data[[1]] %>%
    as_tibble() %>%
    mutate(
      missing = is.na(
        data_list[[simulation_name]]$data[[rel_missingness_index]][[missing_variable]]
      )
    ) %>%
    ggplot(aes(
      x = .data[[missing_variable]],
      y = .data[[comparison_variable]], 
      col = missing,
      alpha = missing
    )) +
    {if (density) geom_density2d()} +
    {if (!density) geom_point()} +
    geom_vline(xintercept = mean_x, alpha = 0.5, linetype = "longdash") +
    geom_hline(yintercept = mean_y, alpha = 0.5, linetype = "longdash") +
    scale_color_manual(values = c("grey", "red")) +
    scale_alpha_manual(
      values = c(alpha_non_missing, 1)
    ) +
    labs(
      x = var2TeX(missing_variable),
      y = var2TeX(comparison_variable),
      title = TeX(paste(
        "Variable with missingness:", var2TeX(missing_variable, tex = FALSE)
      )),
      subtitle = paste(
        "Missingness reason:",
        str_to_upper(data_list[[simulation_name]]$missing_type),
        "\nRelative missingness:",
        relative_missingness
      )) +
    {if (fixed_coords) coord_fixed()} +
    theme_classic() +
    theme(legend.position = "bottom")
  if (!density & marginal) {
    plot <- ggMarginal(plot, type = "density")
  }
  plot
}

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

# simulate_data_list_normal: Generates a list containing simulated data from a multivariate normal distribution.
# Parameters:
#   - seed: Seed for reproducibility of random data generation.
#   - dataset_size: Number of data points to generate (default = 1000).
# Returns a list 'multi_normal_data_list' containing:
#   - 'raw': A nested list storing the simulated multivariate normal data.
#            It contains 'data' as a list with the generated data points,
#            'missing' e.g. x1 if the variable x1 is missing, 'missing_type' e.g.
#           MCAR,MAR or MNAR, and 'relative_missingness' e.g. 0.1, 0.3 or 0.6
simulate_data_list_normal <- function(seed, dataset_size = 1000) {
  multi_normal <- simulate_multi_normal_dataset(
    mu = mu,
    sigma = sigma,
    target_formula = target_formula,
    target_sd = target_sd,
    n = dataset_size,
    seed = seed
  )
  seed <- seed + 1
  set.seed(seed)
  
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
  
  # MCAR for x2 for each relative missingness
  simulation_list <- lapply(relative_missingness, function(rel_missingness) {
    missingness_x2 <- defMiss(varname = "x2", formula = rel_missingness)
    miss_mat <- genMiss(
      multi_normal_data_list$raw$data[[1]],
      missingness_x2,
      idvars = "id"
    )
    genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
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
  
  # MAR for x2 for each relative missingness
  simulation_list <- lapply(relative_missingness, function(rel_missingness) {
    sample_probs <- as.numeric(
      multi_normal_data_list$raw$data[[1]][["x1"]] > 
        mean(multi_normal_data_list$raw$data[[1]][["x1"]])
    ) + 1
    missing_indicator <- sample(
      seq_along(sample_probs),
      size = round(rel_missingness * length(sample_probs)),
      replace = FALSE,
      prob = sample_probs
    )
    missing_indicator <- seq_along(sample_probs) %in% missing_indicator
    missingness_x2 <- defMiss(varname = "x2", formula = "..missing_indicator", logit.link = FALSE)
    miss_mat <- genMiss(
      multi_normal_data_list$raw$data[[1]],
      missingness_x2,
      idvars = "id"
    )
    genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
  })
  
  # add simulated data to data list
  multi_normal_data_list <- add_to_data_list(
    multi_normal_data_list,
    simulation_name = "mar_x2",
    data = simulation_list,
    missing = "x2",
    missing_type = "mar",
    relative_missingness = relative_missingness
  )
  
  # MNAR for x2 for each relative missingness
  simulation_list2 <- lapply(relative_missingness, function(rel_missingness) {
    sample_probs <- as.numeric(
      multi_normal_data_list$raw$data[[1]][["x5"]] > 
        (mean(multi_normal_data_list$raw$data[[1]][["x5"]]))
    ) * 3 + 1
    missing_indicator <- sample(
      seq_along(sample_probs),
      size = round(rel_missingness * length(sample_probs)),
      replace = FALSE,
      prob = sample_probs
    )
    missing_indicator <- seq_along(sample_probs) %in% missing_indicator
    missingness_x2 <- defMiss(varname = "x5", formula = rel_missingness)
    missingness_x2 <- defMiss(
      missingness_x2,
      varname = "x2",
      formula = "..missing_indicator", 
      logit.link = FALSE
    )
    miss_mat <- genMiss(
      multi_normal_data_list$raw$data[[1]],
      missingness_x2,
      idvars = "id"
    )
    genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
    
  })
  
  # add simulated data to data list
  multi_normal_data_list <- add_to_data_list(
    multi_normal_data_list,
    simulation_name = "mnar_x2",
    data = simulation_list2,
    missing = c("x2", "x5"),
    missing_type = c("mnar", "mcar"),
    relative_missingness = relative_missingness
  )
  
  multi_normal_data_list <- create_miss_ind(multi_normal_data_list)
  # get rid of id for the raw data:
  multi_normal_data_list$raw$data[[1]] <- multi_normal_data_list$raw$data[[1]][, -1]
  multi_normal_data_list
}

######################### Mixed Data set ######################### 

# Simulate the raw data
simulate_mixed_dataset <- function(
    n = 1000,
    seed = 2
) {
  mixed_data_def <- defData(
    varname = "x1", dist = "normal", formula = 0, variance = 1
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x2", dist = "beta", formula = 0.75, variance = 3
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x3", dist = "gamma", formula = 2, variance = 2
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x4", dist = "normal", formula = "x3", variance = 4
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x5", dist = "normal", formula = "x4 * 0.5 + x1", variance = 2
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x6", dist = "categorical", 
    formula = "0.2;0.3;0.5",
    variance = "0;1;2"
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x7", dist = "binary", formula = 0.4
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x8", dist = "binary", formula = 0.6
  )
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "x9", dist = "binary", formula = "0.2 + 0.6 * x8"
  )
  
  
  mixed_data_def <- defData(
    mixed_data_def,
    varname = "y",
    formula = "-2 * x1 + 3 * x2 - 0.5 *x3 + 0.5 * x4 + 1.5 * x5 + x6 - 1.5 * x7 + 1.5 * x8 + 2 * x9",
    dist = "normal",
    variance = 5
  )
  set.seed(seed)
  genData(n, mixed_data_def)
}

mixed_df <- simulate_mixed_dataset()

simulate_data_list_mixed <- function(seed, dataset_size = 1000) {
  mixed_df <- simulate_mixed_dataset(seed = seed, n = dataset_size)
  set.seed(seed + 1)
  
  # Store the raw data in the final data list for the multivariate normal 
  # simulations
  mixed_data_list <- list(
    raw = list(
      data = list(mixed_df),
      missing = NA,
      missing_type = NA,
      relative_missingness = 0
    )
  )
  
  # MCAR for each relative missingness
  simulation_list <- lapply(relative_missingness, function(rel_missingness) {
    missingness <- defMiss(varname = "x2", formula = rel_missingness)
    missingness <- defMiss(
      missingness, varname = "x3", formula = rel_missingness
    )
    missingness <- defMiss(
      missingness, varname = "x5", formula = rel_missingness
    )
    missingness <- defMiss(
      missingness, varname = "x6", formula = rel_missingness
    )
    missingness <- defMiss(
      missingness, varname = "x8", formula = rel_missingness
    )
    missingness <- defMiss(
      missingness, varname = "y", formula = rel_missingness
    )
    
    miss_mat <- genMiss(
      mixed_data_list$raw$data[[1]],
      missingness,
      idvars = "id"
    )
    genObs(mixed_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
  })
  
  # add simulated data to data list
  mixed_data_list <- add_to_data_list(
    mixed_data_list,
    simulation_name = "mcar",
    data = simulation_list,
    missing = c("x2", "x3", "x5", "x6", "x8", "x10"),
    missing_type = rep("mcar", 6),
    relative_missingness = relative_missingness
  )
  
  # MAR for each relative missingness
  simulation_list <- lapply(relative_missingness, function(rel_missingness) {
    missingness_indicator_x2 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      multiplier = 3
    )
    missingness <- defMiss(
      varname = "x2",
      formula = "..missingness_indicator_x2",
      logit.link = FALSE
    )
    
    missingness_indicator_x3 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      seed = 3,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x3",
      formula = "..missingness_indicator_x3",
      logit.link = FALSE
    )
    
    missingness_indicator_x5 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      seed = 5,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x5",
      formula = "..missingness_indicator_x5",
      logit.link = FALSE
    )
    
    missingness_indicator_x6 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x7",
      rel_missingness = rel_missingness,
      seed = 6,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x6",
      formula = "..missingness_indicator_x6",
      logit.link = FALSE
    )
    
    missingness_indicator_x8 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x7",
      rel_missingness = rel_missingness,
      seed = 8,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x8",
      formula = "..missingness_indicator_x8",
      logit.link = FALSE
    )
    
    missingness_indicator_from1 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      seed = 11,
      multiplier = 2.2
    )
    missingness_indicator_from4 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x4",
      rel_missingness = rel_missingness,
      seed = 14,
      multiplier = 2.2
    )
    missingness_positions_from1 <- which(missingness_indicator_from1)
    missingness_positions_from1 <- sample(
      missingness_positions_from1,
      size = round(0.8 * length(missingness_positions_from1))
    )
    missingness_positions_from4 <- which(missingness_indicator_from4)
    missingness_positions_from4 <- sample(
      missingness_positions_from4,
      size = round(0.2 * length(missingness_positions_from4))
    )
    
    combined_indicator <- (
      (seq_along(missingness_indicator_from1) %in% missingness_positions_from1) |
        (seq_along(missingness_indicator_from1) %in% missingness_positions_from4)
    )
    missingness <- defMiss(
      missingness,
      varname = "y",
      formula = "..combined_indicator",
      logit.link = FALSE
    )
    
    miss_mat <- genMiss(
      mixed_data_list$raw$data[[1]],
      missingness,
      idvars = "id"
    )
    genObs(mixed_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
  })
  
  # add simulated data to data list
  mixed_data_list <- add_to_data_list(
    mixed_data_list,
    simulation_name = "mar",
    data = simulation_list,
    missing = c("x2", "x3", "x5", "x6", "x8", "y"),
    missing_type = rep("mar", 6),
    relative_missingness = relative_missingness
  )
  
  # MNAR for each relative missingness
  simulation_list <- lapply(relative_missingness, function(rel_missingness) {
    missingness <- defMiss(
      varname = "x1",
      formula = rel_missingness
    )
    
    missingness_indicator_x2 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x2",
      formula = "..missingness_indicator_x2",
      logit.link = FALSE
    )
    
    
    missingness_indicator_x3 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x2",
      rel_missingness = rel_missingness,
      seed = 3,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x3",
      formula = "..missingness_indicator_x3",
      logit.link = FALSE
    )
    
    missingness_indicator_x5fromx1 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      seed = 5,
      multiplier = 3
    )
    missingness_indicator_x5fromx4 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x4",
      rel_missingness = rel_missingness,
      seed = 55,
      multiplier = 3
    )
    
    missingness_positions_x5fromx1 <- which(missingness_indicator_x5fromx1)
    missingness_positions_x5fromx1 <- sample(
      missingness_positions_x5fromx1,
      size = round(0.8 * length(missingness_positions_x5fromx1))
    )
    missingness_positions_x5fromx4 <- which(missingness_indicator_x5fromx4)
    missingness_positions_x5fromx4 <- sample(
      missingness_positions_x5fromx4,
      size = round(0.2 * length(missingness_positions_x5fromx4))
    )
    
    combined_indicator_x5 <- (
      (seq_along(missingness_indicator_x5fromx1) %in% missingness_positions_x5fromx1) |
        (seq_along(missingness_indicator_x5fromx4) %in% missingness_positions_x5fromx4)
    )
    
    missingness <- defMiss(
      missingness,
      varname = "x5",
      formula = "..combined_indicator_x5",
      logit.link = FALSE
    )
    
    missingness <- defMiss(
      missingness,
      varname = "x7",
      formula = rel_missingness
    )
    
    missingness_indicator_x6 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x7",
      rel_missingness = rel_missingness,
      seed = 6,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x6",
      formula = "..missingness_indicator_x6",
      logit.link = FALSE
    )
    
    missingness_indicator_x8 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x7",
      rel_missingness = rel_missingness,
      seed = 8,
      multiplier = 3
    )
    missingness <- defMiss(
      missingness,
      varname = "x8",
      formula = "..missingness_indicator_x8",
      logit.link = FALSE
    )
    
    missingness_indicator_from1 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x1",
      rel_missingness = rel_missingness,
      seed = 11,
      multiplier = 3
    )
    missingness_indicator_from7 <- create_missingness_indicator(
      mixed_data_list$raw$data[[1]],
      from = "x7",
      rel_missingness = rel_missingness,
      seed = 14,
      multiplier = 3
    )
    combined_indicator <- missingness_indicator_from1 | missingness_indicator_from7
    target_missings <- round(rel_missingness * dim(mixed_data_list$raw$data[[1]])[1])
    set_to_false <- sample(which(combined_indicator), sum(combined_indicator) - target_missings)
    combined_indicator[set_to_false] = FALSE
    missingness <- defMiss(
      missingness,
      varname = "y",
      formula = "..combined_indicator",
      logit.link = FALSE
    )
    
    miss_mat <- genMiss(
      mixed_data_list$raw$data[[1]],
      missingness,
      idvars = "id"
    )
    genObs(mixed_data_list$raw$data[[1]], miss_mat, idvars = "id")[, -1]
  })
  
  # add simulated data to data list
  mixed_data_list <- add_to_data_list(
    mixed_data_list,
    simulation_name = "mnar",
    data = simulation_list,
    missing = c("x1", "x2", "x3", "x5", "x6", "x7", "x8", "y"),
    missing_type = c("mcar", "mnar", "mnar", "mnar","mnar", "mcar", "mnar", "mnar"),
    relative_missingness = relative_missingness
  )
  
  
  mixed_data_list <- create_miss_ind(mixed_data_list)
  # get rid of id for the raw data:
  mixed_data_list$raw$data[[1]] <- mixed_data_list$raw$data[[1]][, -1]
  mixed_data_list
}

