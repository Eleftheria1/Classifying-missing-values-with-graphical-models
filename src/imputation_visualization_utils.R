# Function to compare univariate densities of variables
ggplot_compare_density <- function(
    amelia_obj,
    data_list, 
    var = "x2",
    nominal = FALSE,
    expand_limits = NULL
) {
  rel_missingness <- mean(amelia_obj$missMatrix[, var])
  
  plot_observed <- data.frame(
    observed = amelia_obj$imputations[[1]] %>% 
      filter(!amelia_obj$missMatrix[, var]) %>%
      pull(var)
  )
  
  feature_range <- abs(
    max(plot_observed$observed) - min(plot_observed$observed)
  )
  if (is.null(expand_limits) & (feature_range > 2)) {
    expand_limits <- -1 * feature_range * 0.03
  } else {
    expand_limits <- 0
  }
  
  if (nominal) {
    aggregation_fun <- function(x) {
      as.numeric(
        apply(x, 1, function(r) names(which.max(table(r))))
      )
    }
  } else {
    aggregation_fun <- rowMeans
    # only the first imputation:
    # aggregation_fun <- function(x) x[, 1]
  }
  
  
  
  if (nominal) {
    plot_missing <- data.frame(
      missing_means = do.call(
        aggregation_fun,
        list(
          x = vapply(
            amelia_obj$imputations, 
            FUN = function(df) {
              tmp <- df %>%
                pull(var)
              tmp
            },
            FUN.VALUE = numeric(
              length = dim(amelia_obj$missMatrix)[1]
            )
          )
        )
      )
    )
    plot_observed %>%
      group_by(observed) %>%
      summarise(count = n()) %>% 
      ungroup() %>%
      mutate(count = count/sum(count)) %>%
      rename(values = "observed", rel_count = "count") %>%
      mutate(type = "observed") %>%
      bind_rows(
        plot_missing %>%
          group_by(missing_means) %>%
          summarise(count = n()) %>% 
          ungroup() %>%
          mutate(count = count/sum(count)) %>%
          rename(values = "missing_means", rel_count = "count") %>%
          mutate(type = "imputed")
      ) %>%
      bind_rows(
        tibble(true_vals =  data_list$raw$data[[1]][[var]]) %>%
          group_by(true_vals) %>%
          summarise(count = n()) %>% 
          ungroup() %>%
          mutate(count = count/sum(count)) %>%
          rename(values = "true_vals", rel_count = "count") %>%
          mutate(type = "true_vals")
      ) %>%
      mutate(values = factor(values)) %>%
      ggplot() +
      geom_col(
        aes(x = values, y = rel_count, fill = type),
        position = "dodge"
      ) +
      labs(
        x = paste(var, "- Fraction Missing:", round(rel_missingness, 1)),
        y = "Relative Frequency",
        title = paste("Observed, imputed and true frequencies of", var)
      ) +
      scale_fill_manual(
        values = c("observed" = "#33a3f2", "imputed" = "#eb4b3d", "true_vals" = "#ababab"),
        name = "",
        labels = c("Majority Imputations", "Observed Values", "True Values")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  } else {
    plot_missing <- data.frame(
      missing_means_y = do.call(
        aggregation_fun,
        list(
          x = vapply(
            amelia_obj$imputations, 
            FUN = function(df) {
              tmp <- df %>%
                pull(var)
              density(tmp)$y
            },
            FUN.VALUE = numeric(
              length = length(density(amelia_obj$imputations[[1]][[var]])$y)
            )
          )
        )
      ),
      missing_means_x = do.call(
        aggregation_fun,
        list(
          x = vapply(
            amelia_obj$imputations, 
            FUN = function(df) {
              tmp <- df %>%
                pull(var)
              density(tmp)$x
            },
            FUN.VALUE = numeric(
              length = length(density(amelia_obj$imputations[[1]][[var]])$x)
            )
          )
        )
      )
    )
    ggplot() +
      geom_line(
        data = data.frame(
          x = density(data_list$raw$data[[1]][[var]])$x,
          y = density(data_list$raw$data[[1]][[var]])$y
        ),
        aes(x = x, y = y, col = "true_vals", linetype = "true_vals"),
        linewidth = 0.8,
        alpha = 0.7
      ) +
      geom_line(
        data = data.frame(
          x = density(plot_observed$observed)$x,
          y = density(plot_observed$observed)$y
        ),
        aes(x = x, y = y, col = "observed", linetype = "observed"),
        linewidth = 0.8,
        alpha = 1
      ) +
      geom_line(
        data = data.frame(
          x = plot_missing$missing_means_x,
          y = plot_missing$missing_means_y
        ),
        aes(x = x, y = y, col = "imputed", linetype = "imputed"),
        linewidth = 0.8,
        alpha = 1
      ) +
      expand_limits(x = expand_limits) +
      labs(
        x = paste(var, "- Fraction Missing:", round(rel_missingness, 3)),
        y = "Relative Density",
        title = paste("Observed, imputed and true distribution of", var)
      ) +
      scale_color_manual(
        values = c("observed" = "#33a3f2", "imputed" = "#eb4b3d", "true_vals" = "#4b4c4d"),
        name = "",
        labels = c("Averaged Imputations", "Observed Values", "True Values")
      ) +
      scale_linetype_manual(
        values = c("observed" = "solid", "imputed" = "solid", "true_vals" = "longdash"),
        name = "",
        labels = c("Averaged Imputations", "Observed Values", "True Values")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }
}


library(tictoc)
add_all_compare_density_plots <- function(
    imputed_data_list,
    nominal_features = character()
) {
  tic()
  for (exp in names(imputed_data_list)[-1]) {
    for (i in seq_along(imputed_data_list[[exp]]$data)) {
      for (var in imputed_data_list[[exp]]$missing) {
        imputed_data_list[[exp]]$dens_plots[[var]][[i]] <- ggplot_compare_density(
          amelia_obj = imputed_data_list[[exp]]$amelia_obj[[i]],
          data_list = imputed_data_list,
          var = var,
          nominal = ifelse(var %in% nominal_features, TRUE, FALSE)
        )
      }
    }
    cat("\nEperiment", exp, "finished.")
  }
  time <- toc()
  imputed_data_list
}



#################################################################################
# "-2 * x1 + 3 * x2 - 0.5 *x3 + 0.5 * x4 + 1.5 * x5 + x6 - 1.5 * x7 + 1.5 * x8 + 2 * x9"
#visualize_parameters <- function(
#    parameter_result_df,
#    true_params = c(
#      "(Intercept)" = 0,
#      "x1" = -2,
#      "x2" = 3,
#      "x3" = -0.5,
#      "x4" = 0.5,
#      "x5" = 1.5,
#      "x61" = 1,
#      "x62" = 2,
#      "x7" = -1.5,
#      "x8" = 1.5,
#      "x9" = 2
#    )
#) {
#  plot_data <- parameter_result_df %>%
#    bind_rows(
#      tibble(
#        est = true_params,
#        coefficient = names(true_params),
#        lower_bound = true_params,
#        upper_bound = true_params,
#        type = "original",
#        complete_cases = NA
#      )
#    ) %>%
#    mutate(
#      type = case_when(
#        type == "original" ~ "truth",
#        type == "raw" ~ "full population",
#        type == "complete" ~ "complete cases",
#        TRUE ~ "imputed"
#      )
#    ) %>%
#    mutate(type = fct_relevel(
#      factor(type),
#      c("truth", "full population", "complete cases", "imputed")
#    ))
#  plot_data %>%
#    ggplot() +
#    geom_point(
#      aes(x = type, y = est, color = type)
#    ) +
#    geom_errorbar(
#      aes(x = type, ymin = lower_bound, ymax= upper_bound, color = type),
#      data = plot_data %>% filter(type != "truth")
#    ) +
#    facet_wrap(~coefficient, ncol = 4, scales = "free_y") +
#    scale_color_manual(name = "Estimation", values = c(
#      "#4b4c4d",
#      "#4b6da3",
#      "#33a3f2",
#      "#eb4b3d"
#    )) +
#    geom_hline(
#      aes(yintercept = est), 
#      data = plot_data %>% filter(type == "truth"),
#      col = "#4b4c4d",
#      linetype = "dotted"
#    ) +
#    labs(
#      y = "Coefficients", x = "",
#      caption = paste(
#        "#Complete Cases:",
#        plot_data %>% 
#          filter(type == "complete cases") %>%
#          pull(complete_cases) %>%
#          unique()
#      )
#    ) +
#    theme_minimal() +
#    theme(axis.text.x = element_blank())
#  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#}
#



visualize_parameters <- function(
    parameter_result_df
) {
  plot_data <- parameter_result_df %>%
    mutate(
      type = case_when(
        type == "raw" ~ "full population",
        type == "complete" ~ "complete cases",
        TRUE ~ "imputed"
      )
    ) %>%
    mutate(type = fct_relevel(
      factor(type),
      c("full population", "complete cases", "imputed")
    ))
  
  # Create a tibble with the full population values for each coefficient
  full_population_values <- plot_data %>%
    filter(type == "full population") %>%
    select(coefficient, est) %>%
    rename(full_population = est)
  
  plot_data %>%
    ggplot() +
    geom_point(
      aes(x = type, y = est, color = type)
    ) +
    geom_errorbar(
      aes(x = type, ymin = lower_bound, ymax= upper_bound, color = type),
      data = plot_data
    ) +
    facet_wrap(~coefficient, ncol = 4, scales = "free_y") +
    scale_color_manual(name = "Estimation", values = c(
      "#4b6da3",
      "#33a3f2",
      "#eb4b3d"
    )) +
    labs(
      y = "Coefficients", x = "",
      caption = paste(
        "#Complete Cases:",
        plot_data %>% 
          filter(type == "complete cases") %>%
          pull(complete_cases) %>%
          unique()
      )
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    # Add a geom_hline for each coefficient at the full population value
    geom_hline(
      aes(yintercept = full_population),
      data = full_population_values,
      col = "#4b4c4d",
      linetype = "dotted"
    )
}


############################################################################
#knn
# visualize

# selbe diagnostische plots
knn_density_comparison_plots <- function(
    knn_imputed_data_list, col_name, rel_miss, exp = "mnar"
) {
  rel_miss_ind <- which(
    rel_miss == knn_imputed_data_list[[exp]]$relative_missingness
  )
  if (length(rel_miss_ind) == 0) {
    stop("This Relative missingness does not exist.")
  }
  plot_knn_df <- data.frame(
    missing_means_y = do.call(
      rowMeans,
      list(
        x = vapply(
          knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]], 
          FUN = function(df) {
            tmp <- df %>%
              pull(col_name)
            density(tmp)$y
          },
          FUN.VALUE = numeric(
            length = length(
              density(knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]]$imp1[[col_name]])$y
            )
          )
        )
      )
    ),
    missing_means_x = do.call(
      rowMeans,
      list(
        x = vapply(
          knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]], 
          FUN = function(df) {
            tmp <- df %>%
              pull(col_name)
            density(tmp)$x
          },
          FUN.VALUE = numeric(
            length = length(
              density(knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]]$imp1[[col_name]])$y
            )
          )
        )
      )
    )
  )
  
  knn_imputed_data_list[[exp]]$dens_plots[[col_name]][[rel_miss_ind]] +
    geom_line(
      # data = data.frame(
      #   x = density(
      #     knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]][[col_name]]
      #   )$x,
      #   y = density(
      #     knn_imputed_data_list[[exp]]$knn_obj[[rel_miss_ind]][[col_name]]
      #   )$y
      # ),
      data = data.frame(
        x = plot_knn_df$missing_means_x,
        y = plot_knn_df$missing_means_y
      ),
      aes(x = x, y = y), col = "#f0b74d", linetype = "solid",
      linewidth = 0.8,
      alpha = 1,
      inherit.aes = F
    ) +
    labs(subtitle = paste(
      "KNN imputation displayed in ",
      "<span style='color:",
      "#f0b74d",
      "'>**orange**</span>"
    )) +
    theme(plot.subtitle = ggtext::element_markdown(size = 9))
}



visualize_parameters_knn <- function(
    parameter_results,
    classified_data_list,
    experiment = "mnar",
    rel_missingness = 0.3,
    estimated_missingness = TRUE,
    complete_case_threshold = 5
) {
  rel_missingness_index <- which(
    classified_data_list[[experiment]]$relative_missingness == rel_missingness
  )
  parameter_result_df <- parameter_results[[experiment]][[rel_missingness_index]]
  complete_case_note <- ""
  if (min(parameter_result_df$complete_cases) <= complete_case_threshold) {
    parameter_result_df <- parameter_result_df %>%
      filter(type != "complete")
    complete_case_note <- paste("Due to only", min(parameter_result_df$complete_cases),
        "complete cases the corresponding estimates are not displayed.")
    cat(complete_case_note)
  } else {
    complete_case_note <- ""
  }
  
  if (estimated_missingness) {
    treated_mnar <- str_flatten(
      names(
        classified_data_list[[experiment]]$detected_missingness_type[[rel_missingness_index]][
          classified_data_list[[experiment]]$detected_missingness_type[[rel_missingness_index]] == "mnar"
        ]
      ),
      collapse = ","
    )
  } else {
    treated_mnar <- str_flatten(
      names(
        classified_data_list[[experiment]]$missing_type[
          classified_data_list[[experiment]]$missing_type == "mnar"
        ]
      ),
      collapse = ", "
    )
  }
  treated_mnar <- ifelse(treated_mnar == "", "None", treated_mnar)
  
  # display only the knn imputations if it was performed
  fct_levels <- c("full population", "complete cases", "imputed", "KNN imputed")
  fct_cols <- c(
    "#4b6da3",
    "#33a3f2",
    "#eb4b3d",
    "#f0b74d"
  )
  if (complete_case_note != "") {
    fct_levels <- fct_levels[-2]
    fct_cols <- fct_cols[-2]
  }
  plot_data <- parameter_result_df %>%
    mutate(
      type = case_when(
        type == "raw" ~ "full population",
        type == "complete" ~ "complete cases",
        type == "knn_imputed" ~ "KNN imputed",
        TRUE ~ "imputed"
      )
    ) %>%
    mutate(type = fct_relevel(
      factor(type),
      fct_levels
    ))
  
  if (complete_case_note == "") {
    complete_case_note <- paste(
      "#Complete Cases:",
      plot_data %>% 
        filter(type == "complete cases") %>%
        pull(complete_cases) %>%
        unique()
    )
  }
  
  
  plot_data %>%
    ggplot() +
    geom_point(
      aes(x = type, y = est, color = type)
    ) +
    geom_errorbar(
      aes(x = type, ymin = lower_bound, ymax = upper_bound, color = type),
      data = plot_data %>% filter(type != "truth")
    ) +
    facet_wrap(~coefficient, ncol = 4, scales = "free_y") +
    scale_color_manual(name = "Estimation", values = 
      fct_cols
    ) +
    geom_hline(
      aes(yintercept = est), 
      data = plot_data %>% filter(type == "full population"),
      col = "#4b4c4d",
      linetype = "dotted"
    ) +
    labs(
      y = "Coefficients", x = "",
      title = paste("Experiment:", str_to_upper(experiment), "- Relative missingness:", rel_missingness),
      subtitle = paste("Treated as MNAR:", treated_mnar),
      caption = complete_case_note
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


knn_categorical_freq_table <- function(
    data_list,
    var = "x2",
    exp = "mnar",
    rel_miss = 0.1
) {
  rel_missingness_index <- c(
    "0.1" = 1,
    "0.3" = 2,
    "0.6" = 3
  )[as.character(rel_miss)]
  
  # complete cases
  plot_observed <- data.frame(
    observed = data_list[[exp]]$amelia_obj[[rel_missingness_index]]$imputations[[1]] %>% 
      filter(!data_list[[exp]]$amelia_obj[[rel_missingness_index]]$missMatrix[, var]) %>%
      pull(var)
  )
  
  
  aggregation_fun <- function(x) {
    as.numeric(
      apply(x, 1, function(r) names(which.max(table(r))))
    )
  }
  # amelia imputations
  plot_missing <- data.frame(
    missing_means = do.call(
      aggregation_fun,
      list(
        x = vapply(
          data_list[[exp]]$amelia_obj[[rel_missingness_index]]$imputations, 
          FUN = function(df) {
            tmp <- df %>%
              pull(var)
            tmp
          },
          FUN.VALUE = numeric(
            length = dim(data_list[[exp]]$amelia_obj[[rel_missingness_index]]$missMatrix)[1]
          )
        )
      )
    )
  )
  # knn imputations
  plot_missing_knn <- data.frame(
    missing_means = do.call(
      aggregation_fun,
      list(
        x = vapply(
          data_list[[exp]]$knn_obj[[rel_missingness_index]], 
          FUN = function(df) {
            tmp <- df %>%
              pull(var)
            tmp
          },
          FUN.VALUE = numeric(
            length = dim(data_list[[exp]]$amelia_obj[[rel_missingness_index]]$missMatrix)[1]
          )
        )
      )
    )
  )
  tmp <- plot_observed %>%
    group_by(observed) %>%
    summarise(count = n()) %>% 
    ungroup() %>%
    mutate(count = count/sum(count)) %>%
    rename(values = "observed", rel_count = "count") %>%
    mutate(type = "observed") %>%
    bind_rows(
      plot_missing %>%
        group_by(missing_means) %>%
        summarise(count = n()) %>% 
        ungroup() %>%
        mutate(count = count/sum(count)) %>%
        rename(values = "missing_means", rel_count = "count") %>%
        mutate(type = "imputed_amelia")
    ) %>%
    bind_rows(
      plot_missing_knn %>%
        group_by(missing_means) %>%
        summarise(count = n()) %>% 
        ungroup() %>%
        mutate(count = count/sum(count)) %>%
        rename(values = "missing_means", rel_count = "count") %>%
        mutate(type = "imputed_knn")
    ) %>%
    bind_rows(
      tibble(true_vals =  data_list$raw$data[[1]][[var]]) %>%
        group_by(true_vals) %>%
        summarise(count = n()) %>% 
        ungroup() %>%
        mutate(count = count/sum(count)) %>%
        rename(values = "true_vals", rel_count = "count") %>%
        mutate(type = "true_vals")
    ) %>%
    mutate(values = factor(values)) %>%
    pivot_wider(names_from = type, values_from = rel_count, id_cols = values)
  tmp %>%
    mutate(
      across(
        !contains("val"),
        ~ (.x - tmp$true_vals)
      )
    )
}




