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


# MAR for x2 for each relative missingness
simulation_list <- lapply(relative_missingness, function(rel_missingness) {
  set.seed(2)
  sample_probs <- as.numeric(
    multi_normal_data_list$raw$data[[1]][["x1"]] > 
      mean(multi_normal_data_list$raw$data[[1]][["x1"]])) + 1
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
  genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")
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
simulation_list <- lapply(relative_missingness, function(rel_missingness) {
  set.seed(2)
  sample_probs <- as.numeric(
    multi_normal_data_list$raw$data[[1]][["x2"]] > 
      mean(multi_normal_data_list$raw$data[[1]][["x2"]])) + 1
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
  genObs(multi_normal_data_list$raw$data[[1]], miss_mat, idvars = "id")
})

# hist(simulation_list[[1]]$x1[is.na(simulation_list[[1]]$x2)])

# add simulated data to data list
multi_normal_data_list <- add_to_data_list(
  multi_normal_data_list,
  simulation_name = "mnar_x2",
  data = simulation_list,
  missing = "x2",
  missing_type = "mnar",
  relative_missingness = relative_missingness
)

var2TeX <- function(string, tex = TRUE) {
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
    density = FALSE
  ) {
  rel_missingness_index <- which(
    data_list[[simulation_name]]$relative_missingness == relative_missingness
  )
  missing_variable <- data_list[[simulation_name]]$missing
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
  simulation_name = "mar_x2",
  relative_missingness = 0.1,
  alpha_non_missing = 0.3,
  density = F
)


save(
  multi_normal_data_list,
  file = paste0(here::here(), "/data/multi_normal.RData")
)





