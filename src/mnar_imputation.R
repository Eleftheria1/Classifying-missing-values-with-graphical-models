library(kknn)

source("src/imputation_utils.R")
source("src/imputation_visualization_utils.R")
load("data/imputed_mixed_data.RData")
load("data/imputed_norm_mixed_data.RData")
load("data/imputed_normal_data.RData")
load("data/classified_multi_normal_data_list.RData")
load("data/classified_mixed_data_list.RData")

if ("x10" %in% mixed_data_list$mcar$missing) {
  mixed_data_list$mcar$missing[6] <- "y"
}
if ("x10" %in% imputed_mixed_data_list$mcar$missing) {
  imputed_mixed_data_list$mcar$missing[6] <- "y"
}
if ("x10" %in% imputed_norm_mixed_data_list$mcar$missing) {
  imputed_norm_mixed_data_list$mcar$missing[6] <- "y"
}

###############################################################################
# Visualize densities
vizualized_multi_normal_list <- add_all_compare_density_plots(
  imputed_normal_data_list
)
visualized_mixed_data_list <- add_all_compare_density_plots(
  imputed_mixed_data_list,
  nominal_features = c("x6", "x7", "x8", "x9")
)
visualized_norm_mixed_data_list <- add_all_compare_density_plots(
  imputed_norm_mixed_data_list,
  nominal_features = c("x6", "x7", "x8", "x9")
)
# knn imputation for the ground truth mnar columns
# for truth one should also specify the true neighbors if graph_neigbors is used
# this can be done via the feature_vec_truth argument that takes a list of 
# nodes with each of them with their neighbors like 
# list("x1" = c("x2", "x7"), "x2" = ...)
knn_norm_imputed_mixed_data_list <- overall_impute_knn(
  visualized_norm_mixed_data_list,
  mnar_selection = "truth",
  graph_neighbors = FALSE,
  k = 20,
)
# knn imputation for all the columns that were classified to be MNAR by the
# m-graphs
knn_norm_imputed_mixed_data_list_prediction <- overall_impute_knn(
  visualized_norm_mixed_data_list,
  mnar_selection = classified_mixed_data_list,
  k = 20,
  graph_neighbors = T
)

# for the unnormed visualized mixed data replicate the above
knn_imputed_mixed_data_list <- overall_impute_knn(
  visualized_mixed_data_list,
  mnar_selection = classified_mixed_data_list,
  k = 20
)

#scheint für kleine k besser zu funktionieren
# x4 als input funktioniert nicht weil x4 nicht missing ist!
# funktioniert nicht für categorical variablen aber bei categorical lieber 
# kreuztabelle machen weil man bei Balken e nichts sieht
# x1 und x7 sind mcar. imputation dazu ist amelia sieht aber anders aus als rote linie
# weil wir erste imputation nehmen und nicht mitteln

save_knn_comparison_plots <- function(data_list, data_list_label = NULL) {
  data_list_label <- ifelse(
    is.null(data_list_label), "", paste0(data_list_label, "/")
  )
  for (exp in names(data_list)[-1]) {
    for (i in seq_along(data_list[[exp]]$data)) {
      rel_miss <- data_list[[exp]]$relative_missingness[[i]]
      for (var in c("y", paste0("x", c(1:3, 5)))) {
        if (var %in% names(data_list[[exp]]$dens_plots)) {
          ggsave(
            filename = paste0(
              var, "_density_comparison.png"
            ),
            plot = knn_density_comparison_plots(
              knn_imputed_data_list = data_list,
              col_name = var,
              rel_miss = rel_miss,
              exp = exp
            ),
            device = "png",
            path  = paste0(
              "C:/Users/Eleftheria/OneDrive/Desktop/fs5/Masterarbeit/plots/densities/knn/",
              data_list_label,
              exp, "/",
              rel_miss, "/"
            ),
            bg = "#FFFFFF"
          )
        }
      }
    }
  }
}

save_knn_comparison_plots(
  data_list = knn_norm_imputed_mixed_data_list,
  data_list_label = "norm_mixed"
)

knn_density_comparison_plots(
  knn_imputed_data_list = knn_norm_imputed_mixed_data_list,
  col_name = "x2",
  rel_miss = 0.3,
  exp = "mnar"
)





mixed_data_params_knn <- evaluate_parameters_knn(
  imputed_data_list = knn_imputed_mixed_data_list,
  classified_data_list = classified_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"),
  df = 990,
  estimated_missingness = TRUE
)$parameter_results

norm_mixed_data_params_knn <- evaluate_parameters_knn(
  imputed_data_list = knn_norm_imputed_mixed_data_list,
  classified_data_list = classified_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"),
  df = 990,
  estimated_missingness = TRUE
)$parameter_results

norm_mixed_data_params_knn_prediction <- evaluate_parameters_knn(
  imputed_data_list = knn_norm_imputed_mixed_data_list_prediction,
  classified_data_list = classified_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"),
  df = 990,
  estimated_missingness = TRUE
)$parameter_results


visualize_parameters_knn(
  norm_mixed_data_params_knn_prediction,
  classified_data_list = classified_mixed_data_list,
  experiment = "mnar",
  rel_missingness = 0.1,
  estimated_missingness = TRUE,
  complete_case_threshold = 5
)

visualize_parameters_knn(
  mixed_data_params_knn,
  classified_data_list = classified_mixed_data_list,
  experiment = "mnar",
  rel_missingness = 0.6,
  estimated_missingness = TRUE,
  complete_case_threshold = 5
)


