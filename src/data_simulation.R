######################### Simulate data ######################### 
# For each of the datasets to be simulated we follow the recipe below:
#   1 simulate the raw data
#   2 visualize
#   3 fit a regression
#   4 simulate missing value patterns
#   5 visualize again

# Load the required utilies and functions:
source(paste0(here::here(), "/src/data_simulation_utils.R"))

######################### Multivariate Normal ######################### 
# Here we just perform the computations for an exemplary simulation.
# Here we also want to explore whether the simulation acts as expected.
# However one can easily use the below used functions to replicate multiple such
# simulated datasets

multi_normal_data_list <- simulate_data_list_normal(
  seed = 123,
  dataset_size = 1000
)

# Visualize the raw data
GGally::ggpairs(
  as.data.frame(multi_normal_data_list$raw$data),
  columns = 2:ncol(as.data.frame(multi_normal_data_list$raw$data)),
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
  data = as.data.frame(multi_normal_data_list$raw$data)
) %>%
  summary()

plot_single_missing(
  multi_normal_data_list,
  simulation_name = "mnar_x2",
  relative_missingness = 0.1,
  alpha_non_missing = 0.3,
  comparison_variable = "x5",
  density = F
)

# base graph
pcalg::plot(
  pcalg::pc(
    list(C = cor(
      as.data.frame(multi_normal_data_list$raw$data)),
      n = nrow(as.data.frame(multi_normal_data_list$raw$data))
    ),
    p = ncol(as.data.frame(multi_normal_data_list$raw$data)), 
    alpha = 0.01,
    indepTest = pcalg::gaussCItest
  )
)

# save(
#   multi_normal_data_list,
#   file = paste0(here::here(), "/data/multi_normal.RData")
# )

######################### Mixed Data set ######################### 

mixed_data_list <- simulate_data_list_mixed(
  seed = 123,
  dataset_size = 1000
)

# Visualize the raw data
GGally::ggpairs(
  as.data.frame(mixed_data_list$raw$data),
  columns = 1:ncol(as.data.frame(mixed_data_list$raw$data)),
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
  data = as.data.frame(mixed_data_list$raw$data)
) %>%
  summary()

plot_single_missing(
  mixed_data_list,
  simulation_name = "mnar",
  relative_missingness = 0.1,
  alpha_non_missing = 0.1,
  comparison_variable = "x2",
  density = F,
  missingness_index = 3,
  marginal = T,
  fixed_coords = F
)

# base graph
pcalg::plot(
  pcalg::pc(
    list(C = cor(
      as.data.frame(mixed_data_list$raw$data)),
      n = nrow(as.data.frame(mixed_data_list$raw$data))
    ),
    p = ncol(as.data.frame(mixed_data_list$raw$data)), 
    alpha = 0.01,
    indepTest = pcalg::gaussCItest
  )
)

# save(
#   mixed_data_list,
#   file = paste0(here::here(), "/data/mixed_data.RData")
# )
