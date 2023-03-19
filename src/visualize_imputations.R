
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
  
  plot_missing <- data.frame(
    missing_means = do.call(
      aggregation_fun,
      list(
        x = vapply(
          amelia_obj$imputations, 
          FUN = function(df) {
            df %>%
              # filter(amelia_obj$missMatrix[, var]) %>% 
              pull(var)
          },
          FUN.VALUE = numeric(
            length = dim(amelia_obj$missMatrix)[1]
            # length = sum(amelia_obj$missMatrix[, var])
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
        x = density(plot_missing$missing_means)$x,
        y = density(plot_missing$missing_means)$y
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
ggplot_compare_density(imputed_mixed_data_list$mar$amelia_obj[[3]], mixed_data_list, var = "x2")

ggplot_compare_density(imputed_norm_mixed_data_list2$mar$amelia_obj[[3]], mixed_data_list, var = "x2")

add_all_compare_density_plots <- function(
    imputed_data_list,
    nominal_features = character()
) {
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
  }
  imputed_data_list
}

if ("x10" %in% mixed_data_list$mcar$missing) {
  mixed_data_list$mcar$missing[6] <- "y"
}
if ("x10" %in% imputed_mixed_data_list$mcar$missing) {
  imputed_mixed_data_list$mcar$missing[6] <- "y"
}

vizualized_multi_normal_list <- add_all_compare_density_plots(
  imputed_normal_data_list
)
vizualized_mixed_data_list <- add_all_compare_density_plots(
  imputed_mixed_data_list,
  nominal_features = c("x6", "x7", "x8", "x9")
)
vizualized_mixed_data_list$mar$dens_plots$x8[[2]]
vizualized_mixed_data_list$mar$dens_plots$x2[[1]]

#################################################################################
# "-2 * x1 + 3 * x2 - 0.5 *x3 + 0.5 * x4 + 1.5 * x5 + x6 - 1.5 * x7 + 1.5 * x8 + 2 * x9"
visualize_parameters <- function(
    parameter_result_df,
    true_params = c(
      "(Intercept)" = 0,
      "x1" = -2,
      "x2" = 3,
      "x3" = -0.5,
      "x4" = 0.5,
      "x5" = 1.5,
      "x61" = 1,
      "x62" = 2,
      "x7" = -1.5,
      "x8" = 1.5,
      "x9" = 2
    )
) {
  plot_data <- parameter_result_df %>%
    bind_rows(
      tibble(
        est = true_params,
        coefficient = names(true_params),
        lower_bound = true_params,
        upper_bound = true_params,
        type = "original",
        complete_cases = NA
      )
    ) %>%
    mutate(
      type = case_when(
        type == "original" ~ "truth",
        type == "raw" ~ "full population",
        type == "complete" ~ "complete cases",
        TRUE ~ "imputed"
      )
    ) %>%
    mutate(type = fct_relevel(
      factor(type),
      c("truth", "full population", "complete cases", "imputed")
    ))
  plot_data %>%
    ggplot() +
    geom_point(
      aes(x = type, y = est, color = type)
    ) +
    geom_errorbar(
      aes(x = type, ymin = lower_bound, ymax= upper_bound, color = type),
      data = plot_data %>% filter(type != "truth")
    ) +
    facet_wrap(~coefficient, ncol = 4, scales = "free_y") +
    scale_color_manual(name = "Estimation", values = c(
      "#4b4c4d",
               "#4b6da3",
               "#33a3f2",
               "#eb4b3d"
    )) +
    geom_hline(
      aes(yintercept = est), 
      data = plot_data %>% filter(type == "truth"),
      col = "#4b4c4d",
      linetype = "dotted"
    ) +
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
    theme(axis.text.x = element_blank())
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

visualize_parameters(mixed_data_params$mcar[[2]])