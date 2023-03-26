library(Amelia)
library(clarify)
library(tidyverse)
library(mice)
library(bestNormalize)

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

#save(imputed_mixed_data_list, file = "data/imputed_mixed_data.RData")

imputed_normal_data_list <- overall_impute(
  data_list = multi_normal_data_list
)

factorize_columns <- function(df, columns = c("x6")) {
  for (col in columns) {
    if (length(unique(df[[col]])) > 1) {
      df[[col]] <- factor(df[[col]])
    }
  }
  df
}

evaluate_parameters <- function(
    imputed_data_list,
    form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"),
    df = 990
) {
  model_list <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  pooled_results <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  parameter_result <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  
  # raw data --------------------------------------------
  model_raw <- lm(
    formula = form,
    data = factorize_columns(imputed_data_list$raw$data[[1]])
  )
  print(summary(model_raw))
  raw_results <- model_raw$coefficients %>%
    as_tibble() %>%
    select("est" = "value") %>%
    mutate(
      coefficient = rownames(confint(model_raw, level = 0.95)),
      lower_bound = confint(model_raw, level = 0.95)[, 1],
      upper_bound = confint(model_raw, level = 0.95)[, 2],
      type = "raw"
    )
  
  #simulated_coeffs <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  for (exp in names(imputed_data_list)[-1]) {
    for (i in seq_along(imputed_data_list[[exp]]$amelia_obj)) {
      # imputed data --------------------------------------------
      model_list[[exp]][[i]] <- with.amelia(
        factorize_columns(imputed_data_list[[exp]]$amelia_obj[[i]]),
        lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
      )
      pooled_results[[exp]][[i]] <- summary(pool(as.mira(model_list[[exp]][[i]])))
        
      matrix <- matrix(
        c(pooled_results[[exp]][[i]][,2]
        - qt(p = 0.975, df = df) * pooled_results[[exp]][[i]][,3],
        pooled_results[[exp]][[i]][,2]
        + qt(p = 0.975, df = df) * pooled_results[[exp]][[i]][,3]),
        ncol = 2
      )
      
      #pooled_results[[exp]][[i]] <- list(
      #  pooled = pooled_results[[exp]][[i]]$pooled,
      #  glanced = pooled_results[[exp]][[i]]$glanced
      #)
      pooled_results[[exp]][[i]] <- bind_cols(pooled_results[[exp]][[i]], matrix)
      colnames(pooled_results[[exp]][[i]])[7] <- "lower_bound"
      colnames(pooled_results[[exp]][[i]])[8] <- "upper_bound"
      
      # complete case data --------------------------------------------
      model_complete <- lm(
        formula = form,
        data = factorize_columns(na.omit(imputed_data_list[[exp]]$data[[i]]))
      )
      complete_results <- model_complete$coefficients %>%
        as_tibble() %>%
        select("est" = "value") %>%
        mutate(
          coefficient = rownames(confint(model_complete, level = 0.95)),
          lower_bound = confint(model_complete, level = 0.95)[, 1],
          upper_bound = confint(model_complete, level = 0.95)[, 2],
          type = "complete"
        )
      
      # Combine Results
      complete_cases <- nrow(na.omit(imputed_data_list[[exp]]$data[[i]]))
      
      parameter_result[[exp]][[i]] <- pooled_results[[exp]][[i]] %>%
        as_tibble() %>%
        select(
          est = "estimate",
          coefficient = "term",
          lower_bound,
          upper_bound
        ) %>%
        mutate(type = "imputed") %>%
        bind_rows(complete_results) %>%
        bind_rows(raw_results) %>%
        mutate(
          complete_cases = complete_cases,
          coefficient = str_replace_all(
            coefficient, "factor\\(([a-z10-9]+)\\)([0-9])", "\\1\\2"
          )
        ) 
    }
  }
  list(
    model_list = model_list,
    pooled_results = pooled_results,
    parameter_results = parameter_result
  )
}



mixed_data_params <- evaluate_parameters(
  imputed_data_list = imputed_mixed_data_list,
  form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9")
)$parameter_results






# Optional Normalization -------------------------------------------------

# sieht aus als ob die normalisierungs nichts bringt weder
# bei kerneldensity noch bei koeffizienten


normalize_1col <- function(df, col = "x2") {
  bn_obj <- bestNormalize(df[[col]], warn = TRUE, standardize = FALSE)
  df[[col]] <- bn_obj$x.t
  list(bn_obj, df)
}

normalize_columns <- function(data_list, columns = c("x2")) {
  for (exp in names(data_list)[-1]) {
    for (i in seq_along(data_list[[exp]]$data)) {
      data_list[[exp]]$bn_objects[[i]] = list()
      for (col_index in seq_along(columns)) {
        tmp_list <- normalize_1col(
          data_list[[exp]]$data[[i]], col = columns[col_index]
        )
        data_list[[exp]]$bn_objects[[i]] <- append(
          data_list[[exp]]$bn_objects[[i]], list(tmp_list[[1]])
        )
        names(data_list[[exp]]$bn_objects[[i]])[col_index] <- columns[col_index]
        data_list[[exp]]$data[[i]] <- tmp_list[[2]]
      }
    }
  }
  data_list
}

# normalized_mixed_data_list <- normalize_columns(mixed_data_list)

# imputed_norm_mixed_data_list <- overall_impute(
#   data_list = normalized_mixed_data_list,
#   nominals = c("x6", "x7", "x8", "x9"),
#   sqrts = c("x3")
# ) 

retransform_1col <- function(df, bn_obj, coln="x2") {
  df[[coln]] <- bestNormalize:::predict.bestNormalize(
    bn_obj[[coln]], newdata = df[[coln]], inverse = TRUE
  )
  df
}

retransform_columns <- function(imputed_data_list, columns = c("x2")) {
  for (exp in names(imputed_data_list)[-1]) {
    for (i in seq_along(imputed_data_list[[exp]]$data)) {
      for (imp_index in seq_along(imputed_data_list[[exp]]$amelia_obj[[i]]$imputations)) {
        for (col_index in seq_along(columns)) {
          imputed_data_list[[exp]]$amelia_obj[[i]]$imputations[[imp_index]] <- retransform_1col(
            imputed_data_list[[exp]]$amelia_obj[[i]]$imputations[[imp_index]],
            imputed_data_list[[exp]]$bn_objects[[i]][columns[col_index]],
            coln = columns[col_index]
          )
        }
      }
    }
  }
  imputed_data_list
}

# imputed_norm_mixed_data_list <- retransform_columns(imputed_norm_mixed_data_list)

