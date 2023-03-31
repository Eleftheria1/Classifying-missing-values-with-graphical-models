# Evaluate the classification performance

# Load the required utilies and functions:
library(tidyverse)
source(paste0(here::here(), "/src/data_simulation_utils.R"))
source("src/graphical_model_utils.R")
library(progress)
set.seed(123)

evaluate_simulated_classification <- function(
  data_list_type = c("normal", "mixed"),
  dataset_size = 1000,
  replications = 100,
  type = c("pc", "mvpc")
) {
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent)  Time remaining: :eta",
    total = replications
  )
  data_list_type <- match.arg(data_list_type)
  type <- match.arg(type)
  simulation_function <- ifelse(
    data_list_type == "normal",
    simulate_data_list_normal,
    simulate_data_list_mixed
  )
  result_df <- tibble(
    "experiment" = character(),
    "variable" = character(),
    "missingness_type" = character(),
    "relative_missingness" = numeric(),
    "classification" = character(),
    "replication" = numeric()
  )
  pb$tick(0)
  for (i in seq(replications)) {
    pb$tick(1)
    data_list <- do.call(
      simulation_function,
      args = list("seed" = i, "dataset_size" = dataset_size)
    )
    data_list <- fit_and_visualize_graphs(
      data_list = data_list,
      type = type,
      plot_graph = FALSE
    )
    data_list <- add_classification(data_list)
    for (experiment in names(data_list)[2:length(data_list)]) {
      for (rel_missing_i in seq_along(data_list[[experiment]]$relative_missingness)) {
        for (variable_i in seq_along(data_list[[experiment]]$missing)) {
          row_addition <- tibble(
            "experiment" = experiment,
            "variable" = data_list[[experiment]]$missing[[variable_i]],
            "missingness_type" = data_list[[experiment]]$missing_type[[variable_i]],
            "relative_missingness" = data_list[[experiment]]$relative_missingness[[rel_missing_i]],
            "classification" = data_list[[experiment]]$detected_missingness_type[[rel_missing_i]][[variable_i]],
            "replication" = i
          )
          result_df <- bind_rows(result_df, row_addition)
        }
      }
    }
  }
  result_df %>% mutate("correctly_classified" = (classification == missingness_type))
}

test <- evaluate_simulated_classification(
  data_list_type = "mixed",
  dataset_size = 500,
  type = "mvpc"
)

test %>%
  mutate(
    var_type = paste0(variable, " ", str_to_upper(missingness_type)),
    relative_missingness = as.factor(relative_missingness),
    experiment = str_to_upper(str_replace(experiment, "_", " "))
  ) %>%
  group_by(experiment, var_type, relative_missingness) %>%
  summarise(
    correctly_classified = mean(correctly_classified),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = relative_missingness,
      y = var_type,
      fill = correctly_classified, 
      label = correctly_classified
    )
  ) +
  geom_tile(col = "white") +
  geom_text(col = "white", size = 5) +
  labs(x = "Relative Missingness", y = "", fill = "Correctly Classified") +
  facet_wrap(~experiment, scales = "free") +
  theme_minimal()

save(
   test,
   file = paste0(here::here(), "/data/mvpc_results_mixed_500n.RData")
 )


