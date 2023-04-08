library(kknn)

source("src/imputation_utils.R")
source("src/imputation_visualization_utils.R")
load("data/visualized_mixed_data.RData")
load("data/visualized_norm_mixed_data.RData")
load("data/classified_multi_normal_data_list.RData")
load("data/classified_mixed_data_list.RData")

# knn imputation for the ground truth mnar columns
knn_norm_imputed_mixed_data_list <- overall_impute_knn(
  visualized_norm_mixed_data_list,
  mnar_selection = "truth",
  k = 3,
)
# knn imputation for all the columns that were classified to be MNAR by the
# m-graphs
knn_norm_imputed_mixed_data_list <- overall_impute_knn(
  visualized_norm_mixed_data_list,
  mnar_selection = classified_mixed_data_list,
  k = 3,
)

# for the unnormed visualized mixed data replicate the above
knn_imputed_mixed_data_list <- overall_impute_knn(
  visualized_mixed_data_list,
  mnar_selection = classified_mixed_data_list,
  k = 3
)

#scheint für kleine k besser zu funktionieren
# x4 als input funktioniert nicht weil x4 nicht missing ist!
# funktioniert nicht für categorical variablen aber bei categorical lieber 
# kreuztabelle machen weil man bei Balken e nichts sieht
# x1 und x7 sind mcar. imputation dazu ist amelia sieht aber anders aus als rote linie
# weil wir erste imputation nehmen und nicht mitteln

knn_density_comparison_plots(
  knn_imputed_data_list = knn_norm_imputed_mixed_data_list,
  col_name = "x2",
  rel_miss = 0.6,
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


visualize_parameters_knn(
  norm_mixed_data_params_knn,
  classified_data_list = classified_mixed_data_list,
  experiment = "mnar",
  rel_missingness = 0.3,
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


