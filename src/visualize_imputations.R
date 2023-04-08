load("data/imputed_mixed_data.RData")
load("data/imputed_norm_mixed_data.RData")
load("data/imputed_normal_data.RData")
source("src/imputation_utils.R")
source("src/imputation_visualization_utils.R")

if ("x10" %in% mixed_data_list$mcar$missing) {
  mixed_data_list$mcar$missing[6] <- "y"
}
if ("x10" %in% imputed_mixed_data_list$mcar$missing) {
  imputed_mixed_data_list$mcar$missing[6] <- "y"
}
if ("x10" %in% imputed_norm_mixed_data_list$mcar$missing) {
  imputed_norm_mixed_data_list$mcar$missing[6] <- "y"
}
##################################################################
# Visualize coefficients
mixed_data_params <- evaluate_parameters(
  imputed_data_list = imputed_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9")
)$parameter_results

norm_mixed_data_params <- evaluate_parameters(
  imputed_data_list = imputed_norm_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9")
)$parameter_results


visualize_parameters(norm_mixed_data_params$mnar[[1]])
visualize_parameters(mixed_data_params$mnar[[1]])
visualize_parameters(mixed_data_params$mnar[[2]])
visualize_parameters(mixed_data_params$mnar[[3]])
visualize_parameters(mixed_data_params$mar[[1]])
visualize_parameters(mixed_data_params$mar[[2]])
visualize_parameters(mixed_data_params$mar[[3]])
visualize_parameters(mixed_data_params$mcar[[1]])
visualize_parameters(mixed_data_params$mcar[[2]])
visualize_parameters(mixed_data_params$mcar[[3]])

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

#save(visualized_mixed_data_list, file = "data/visualized_mixed_data.RData")
#save(visualized_norm_mixed_data_list, file = "data/visualized_norm_mixed_data.RData")

ggplot_compare_density(imputed_mixed_data_list$mnar$amelia_obj[[3]], mixed_data_list, var = "x2")
ggplot_compare_density(imputed_norm_mixed_data_list$mnar$amelia_obj[[3]], mixed_data_list, var = "x2", nominal = F)


# Variable x_3
visualized_mixed_data_list$mnar$dens_plots$x1[[1]]
visualized_mixed_data_list$mnar$dens_plots$x1[[2]]
visualized_mixed_data_list$mnar$dens_plots$x1[[3]]

# Variable x_7
visualized_mixed_data_list$mnar$dens_plots$x7[[1]]
visualized_mixed_data_list$mnar$dens_plots$x7[[2]]
visualized_mixed_data_list$mnar$dens_plots$x7[[3]]

# Variable x_2
visualized_norm_mixed_data_list$mcar$dens_plots$x2[[1]]
visualized_norm_mixed_data_list$mcar$dens_plots$x2[[2]]
visualized_norm_mixed_data_list$mcar$dens_plots$x2[[3]]
visualized_norm_mixed_data_list$mar$dens_plots$x2[[1]]
visualized_norm_mixed_data_list$mar$dens_plots$x2[[2]]
visualized_norm_mixed_data_list$mar$dens_plots$x2[[3]]
visualized_norm_mixed_data_list$mnar$dens_plots$x2[[1]]
visualized_norm_mixed_data_list$mnar$dens_plots$x2[[2]]
visualized_norm_mixed_data_list$mnar$dens_plots$x2[[3]]

visualized_mixed_data_list$mcar$dens_plots$x2[[1]]
visualized_mixed_data_list$mcar$dens_plots$x2[[2]]
visualized_mixed_data_list$mcar$dens_plots$x2[[3]]
visualized_mixed_data_list$mar$dens_plots$x2[[1]]
visualized_mixed_data_list$mar$dens_plots$x2[[2]]
visualized_mixed_data_list$mar$dens_plots$x2[[3]]
visualized_mixed_data_list$mnar$dens_plots$x2[[1]]
visualized_mixed_data_list$mnar$dens_plots$x2[[2]]
visualized_mixed_data_list$mnar$dens_plots$x2[[3]]

# Variable x_3
visualized_mixed_data_list$mcar$dens_plots$x3[[1]]
visualized_mixed_data_list$mcar$dens_plots$x3[[2]]
visualized_mixed_data_list$mcar$dens_plots$x3[[3]]
visualized_mixed_data_list$mar$dens_plots$x3[[1]]
visualized_mixed_data_list$mar$dens_plots$x3[[2]]
visualized_mixed_data_list$mar$dens_plots$x3[[3]]
visualized_mixed_data_list$mnar$dens_plots$x3[[1]]
visualized_mixed_data_list$mnar$dens_plots$x3[[2]]
visualized_mixed_data_list$mnar$dens_plots$x3[[3]]


# Variable x_5
visualized_mixed_data_list$mcar$dens_plots$x5[[1]]
visualized_mixed_data_list$mcar$dens_plots$x5[[2]]
visualized_mixed_data_list$mcar$dens_plots$x5[[3]]
visualized_mixed_data_list$mar$dens_plots$x5[[1]]
visualized_mixed_data_list$mar$dens_plots$x5[[2]]
visualized_mixed_data_list$mar$dens_plots$x5[[3]]
visualized_mixed_data_list$mnar$dens_plots$x5[[1]]
visualized_mixed_data_list$mnar$dens_plots$x5[[2]]
visualized_mixed_data_list$mnar$dens_plots$x5[[3]]

# Variable x_6
visualized_mixed_data_list$mcar$dens_plots$x6[[1]]
visualized_mixed_data_list$mcar$dens_plots$x6[[2]]
visualized_mixed_data_list$mcar$dens_plots$x6[[3]]
visualized_mixed_data_list$mar$dens_plots$x6[[1]]
visualized_mixed_data_list$mar$dens_plots$x6[[2]]
visualized_mixed_data_list$mar$dens_plots$x6[[3]]
visualized_mixed_data_list$mnar$dens_plots$x6[[1]]
visualized_mixed_data_list$mnar$dens_plots$x6[[2]]
visualized_mixed_data_list$mnar$dens_plots$x6[[3]]

visualized_norm_mixed_data_list$mcar$dens_plots$x6[[1]]
visualized_norm_mixed_data_list$mcar$dens_plots$x6[[2]]
visualized_norm_mixed_data_list$mcar$dens_plots$x6[[3]]
visualized_norm_mixed_data_list$mar$dens_plots$x6[[1]]
visualized_norm_mixed_data_list$mar$dens_plots$x6[[2]]
visualized_norm_mixed_data_list$mar$dens_plots$x6[[3]]
visualized_norm_mixed_data_list$mnar$dens_plots$x6[[1]]
visualized_norm_mixed_data_list$mnar$dens_plots$x6[[2]]
visualized_norm_mixed_data_list$mnar$dens_plots$x6[[3]]

# Variable x_8
visualized_mixed_data_list$mcar$dens_plots$x8[[1]]
visualized_mixed_data_list$mcar$dens_plots$x8[[2]]
visualized_mixed_data_list$mcar$dens_plots$x8[[3]]
visualized_mixed_data_list$mar$dens_plots$x8[[1]]
visualized_mixed_data_list$mar$dens_plots$x8[[2]]
visualized_mixed_data_list$mar$dens_plots$x8[[3]]
visualized_mixed_data_list$mnar$dens_plots$x8[[1]]
visualized_mixed_data_list$mnar$dens_plots$x8[[2]]
visualized_mixed_data_list$mnar$dens_plots$x8[[3]]


