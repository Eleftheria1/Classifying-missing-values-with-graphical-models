######################### Simulate data ######################### 
# For each of the datasets to be simulated we follow the recipe below:
#   1 simulate the raw data
#   2 visualize
#   3 fit a regression
#   4 add missing values
#   5 visualize again

# Load the packages required:
library(tidyverse)
library(ggExtra)
library(latex2exp)
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
  mu = c(0, 5, 10, -1, 3),
  sigma = c(1, 2, 3, 0.3, 1),
  target_formula = "-3 * x1 + 5 * x2 - x3 + 5 * x4 + 2 * x5",
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
  set.seed(2)
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

# hist(simulation_list[[1]]$x1[is.na(simulation_list[[1]]$x2)])

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
  set.seed(2)
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

 #hist(simulation_list[[1]]$x1[is.na(simulation_list[[1]]$x2)])

# add simulated data to data list
# multi_normal_data_list <- add_to_data_list(
#   multi_normal_data_list,
#   simulation_name = "mnar_x2",
#   data = simulation_list,
#   missing = "x2",
#   missing_type = "mnar",
#   relative_missingness = relative_missingness
# )

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

plot_single_missing <- function(
    data_list,
    simulation_name, 
    relative_missingness,
    comparison_variable = "x1",
    alpha_non_missing = 1,
    density = FALSE,
    missingness_index = 1
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
    coord_fixed() +
    theme_classic() +
    theme(legend.position = "bottom")
  if (!density) {
    plot <- ggMarginal(plot, type = "density")
  }
  plot
}

plot_single_missing(
  multi_normal_data_list,
  simulation_name = "mnar_x2",
  relative_missingness = 0.1,
  alpha_non_missing = 0.3,
  comparison_variable = "x5",
  density = F
)

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

multi_normal_data_list <- create_miss_ind(multi_normal_data_list)
# get rid of id for the raw data:
multi_normal_data_list$raw$data[[1]] <- multi_normal_data_list$raw$data[[1]][, -1]
# multi_normal_data_list$raw$data[[1]]


# save(
#   multi_normal_data_list,
#   file = paste0(here::here(), "/data/multi_normal.RData")
# )


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
    varname = "x2", dist = "beta", formula = 0.5, variance = 2
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
# Visualize the raw data
GGally::ggpairs(
  mixed_df,
  columns = 2:ncol(mixed_df),
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
  data = mixed_df %>% select(-id)
) %>%
  summary()

# base graph
pcalg::plot(
  pcalg::pc(
    list(C = cor(mixed_df[, 2:11]), n =1000),
    p = 10, alpha = 0.01,
    indepTest = pcalg::gaussCItest
  )
)

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
  
  set.seed(2)
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

# MAR for each relative missingness
simulation_list <- lapply(relative_missingness, function(rel_missingness) {
  missingness_indicator_x2 <- create_missingness_indicator(
    mixed_data_list$raw$data[[1]],
    from = "x1",
    rel_missingness = rel_missingness
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
    seed = 3
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
    seed = 5
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
    seed = 6
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
    seed = 8
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
  )
  missingness_indicator_from4 <- create_missingness_indicator(
    mixed_data_list$raw$data[[1]],
    from = "x4",
    rel_missingness = rel_missingness,
    seed = 14
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

  set.seed(2)
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

plot_single_missing(
  mixed_data_list,
  simulation_name = "mar",
  relative_missingness = 0.1,
  alpha_non_missing = 0.3,
  comparison_variable = "x1",
  density = F,
  missingness_index = 6
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
    rel_missingness = rel_missingness
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
    seed = 3
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
    seed = 5
  )
  missingness_indicator_x5fromx4 <- create_missingness_indicator(
    mixed_data_list$raw$data[[1]],
    from = "x4",
    rel_missingness = rel_missingness,
    seed = 55
  )
  combined_indicator_x5 <- missingness_indicator_x5fromx1 | missingness_indicator_x5fromx4
  target_missings <- round(rel_missingness * dim(mixed_data_list$raw$data[[1]])[1])
  set_to_false <- sample(which(combined_indicator_x5), sum(combined_indicator_x5) - target_missings)
  combined_indicator_x5[set_to_false] = FALSE
  
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
    seed = 6
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
    seed = 8
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
    seed = 11
  )
  missingness_indicator_from7 <- create_missingness_indicator(
    mixed_data_list$raw$data[[1]],
    from = "x7",
    rel_missingness = rel_missingness,
    seed = 14
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
  
  set.seed(2)
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

### Some debugging code
# View(round(cor(mixed_data_list$mar$data[[2]], use = "pairwise.complete.obs"), 2))
filled_cor <- cor(mixed_data_list$mnar$data[[2]], use = "pairwise.complete.obs")
filled_cor[is.na(filled_cor)] <- 0
pcalg::plot(
  pcalg::pc(
    list(
      C = filled_cor,
      n = 1000
    ),
    p = ncol(mixed_data_list$mnar$data[[2]]), alpha = 0.01,
    indepTest = pcalg::gaussCItest
  )
)

# save(
#   mixed_data_list,
#   file = paste0(here::here(), "/data/mixed_data.RData")
# )
