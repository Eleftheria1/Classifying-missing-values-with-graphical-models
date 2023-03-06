library(Amelia)
library(clarify)
library(tidyverse)
library(mice)

base_datasets <- c("mixed_data.RData", "multi_normal.RData")
all_data_files <- list.files(path = "data")
experiment_results <- all_data_files[!(all_data_files %in% base_datasets)]

load(paste0("data/", base_datasets[1]))
load(paste0("data/", base_datasets[2]))

overall_impute <- function(data_list, nominals = NULL, m = 5, p2s = 0, sqrts = NULL) {
  for (exp in names(data_list)[-1]) {
    for (i in seq_along(data_list[[exp]]$data)) {
      data_list[[exp]]$amelia_obj[[i]] <- amelia(
        data_list[[exp]]$data[[i]] %>%
          select(colnames(data_list$raw$data[[1]])), 
        m = m, 
        noms = nominals,
        p2s = p2s,
        sqrts = sqrts
      )
    }
  }
  data_list
}

set.seed(123)
imputed_mixed_data_list <- overall_impute(
  data_list = mixed_data_list,
  nominals = c("x6", "x7", "x8", "x9"),
  sqrts = c("x3")
)

imputed_normal_data_list <- overall_impute(
  data_list = multi_normal_data_list
)

model_list <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
pooled_results <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
confidence_intervals <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))

#simulated_coeffs <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
for (exp in names(imputed_mixed_data_list)[-1]) {
   for (i in seq_along(imputed_mixed_data_list[[exp]]$amelia_obj)){
         model_list[[exp]][[i]] <- with.amelia(
         imputed_mixed_data_list[[exp]]$amelia_obj[[i]],
         lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9))
         #simulated_coeffs[[exp]][[i]] <- misim(model_list[[exp]][[i]], n = 200)
         pooled_results[[exp]][[i]] <- pool(as.mira(model_list[[exp]][[i]]))
         matrix <- matrix(c(pooled_results[[exp]][[i]]$pooled$estimate
         - qt(p = 0.975, df = 990) * pooled_results[[exp]][[i]]$pooled$t,
         pooled_results[[exp]][[i]]$pooled$estimate
         - qt(p = 0.975, df = 990) * pooled_results[[exp]][[i]]$pooled$t), ncol = 2)
         pooled_results[[exp]][[i]] <- bind_cols(pooled_results[[exp]][[i]]$pooled, matrix)
         colnames(pooled_results[[exp]][[i]])[12] <- "lower_bound"
         colnames(pooled_results[[exp]][[i]])[13] <- "upper_bound"
  }
}

#additional information such as aic or p values
pooled_results$mcar[[3]]$glanced #ist weg wegen konfidenzintervallen
#pooled estimates and variances
pooled_results$mcar[[3]]

#confidence intervals 
t_score <- qt(p = 0.975, df = 990)
conf_lower <- pooled_results$mcar[[3]]$pooled$estimate - t_score * pooled_results$mcar[[3]]$pooled$t
conf_upper <- pooled_results$mcar[[3]]$pooled$estimate + t_score * pooled_results$mcar[[3]]$pooled$t






#transformation notwendig sonst imputation scheiße
#sqrt funktioniert ganz gut bei x3 aber nicht bei x2

# Function to compare univariate densities of variables
# !TODO treat nominal variables with majority class
ggplot_compare_density <- function(
    amelia_obj,
    data_list, 
    var = "x2",
    nominal = FALSE,
    expand_limits = NULL # wenn sd größer 1 dann -0.5 bei numerisch
) {
  rel_missingness <- mean(amelia_obj$missMatrix[, var])
  
  plot_observed <- data.frame(
    observed = amelia_obj$imputations[[1]] %>% 
      filter(!amelia_obj$missMatrix[, var]) %>%
      pull(var)
  )
  
  plot_missing <- data.frame(
    missing_means = rowMeans(
      vapply(
        amelia_obj$imputations, 
        FUN = function(df) {
          df %>% filter(amelia_obj$missMatrix[, var]) %>% pull(var)
        },
        FUN.VALUE = numeric(
          length = sum(amelia_obj$missMatrix[, var])
        )
      )
    )
  )
  
  ggplot() +
    geom_density(
      data = data_list$raw$data[[1]],
      aes(x = .data[[var]], col = "true_vals", linetype = "true_vals"), 
      linewidth = 0.8,
      alpha = 0.7
    ) +
    geom_density(
      data = plot_observed,
      aes(x = observed, col = "observed", linetype = "observed"), 
      linewidth = 0.8,
      alpha = 1
    ) +
    geom_density(
      data = plot_missing,
      aes(x = missing_means, col = "imputed", linetype = "imputed"),
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
ggplot_compare_density(imputed_mixed_data_list$mar$amelia_obj[[3]], mixed_data_list, var = "x3")
ggplot_compare_density(a.out, multi_normal_data_list, var = "x2")

add_all_compare_density_plots <- function(imputed_data_list) {
  for (exp in names(imputed_data_list)[-1]) {
    for (i in seq_along(imputed_data_list[[exp]]$data)) {
      for (var in imputed_data_list[[exp]]$missing) {
        imputed_data_list[[exp]]$dens_plots[[var]][[i]] <- ggplot_compare_density(
          amelia_obj = imputed_data_list[[exp]]$amelia_obj[[i]],
          data_list = imputed_data_list,
          var = var,
          expand_limits = -0.5
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

vizualized_multi_normal_list <- add_all_compare_density_plots(imputed_normal_data_list)
vizualized_mixed_data_list <- add_all_compare_density_plots(imputed_mixed_data_list)


# -------------------------------------------

#transformation notwendig sonst imputation scheiße
#sqrt funktioniert ganz gut bei x3 aber nicht bei x2
plot(imputed_mixed_data_list$mar$amelia_obj[[3]])
lines(density(mixed_data_list$raw$data[[1]]$x3))
legend("topright", legend = c("Mean Imputations", 
                              "Observed Values"))


plot(imputed_mixed_data_list$mar$amelia_obj[[3]], which.vars = 2)
lines(density(mixed_data_list$raw$data[[1]]$x2))

disperse(imputed_mixed_data_list$mar$amelia_obj[[3]], dims = 2, m = 5)











data <- multi_normal_data_list$mcar_x2$data[[3]]
a.out <- amelia(data, m = 5, noms = "missing_x2", p2s = 2)

# Fit model to each dataset
model.list <- with(a.out, lm(y ~ x1 + x2 + x3 + x4 +x5))
# Simulate coefficients
si <- misim(model.list, n = 200)
si$fit
# compute the average marginal effect of x2: 4.888511
effect <- sim_ame(si, var = "x2")
effect
#compare observed effect of x2 4.9305
lm(formula = y ~ . - missing_x2, data = data)
# and true value of x2 = 4.9344 
lm(formula = y ~ . , data = multi_normal_data_list$raw$data[[1]])

data.frame(x2_coef = si$sim.coefs[,3]) %>%
  mutate(imputation = si$imp) %>%
  group_by(imputation) %>%
  summarize(group_mean = mean(x2_coef)) %>%
  pull(group_mean) %>%
  mean()


#red means mean imputations blue means observed values
# lines() funktioniert so nur wenn es nur eine missing variable gibt

#plot(a.out, which.vars = 1:7, lwd = 2)



compare.density(a.out, var = "x2", lwd = 2)
lines(density(multi_normal_data_list$raw$data[[1]]$x2), lwd = 2, lty = 2, col = "black")

test <- overimpute(a.out, var = "x5")
test$lower.overimputed
# check convergence and starting points
disperse(a.out, dims = 2, m = 5)


