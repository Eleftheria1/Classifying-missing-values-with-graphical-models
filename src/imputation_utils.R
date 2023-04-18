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
      imputation_length <- length(
        imputed_data_list[[exp]]$amelia_obj[[i]]$imputations
      )
      tmp_amelia_obj_list <- lapply((1:(imputation_length / 5)) - 1, function(j) {
        list(imputations = lapply(seq(j * 5 + 1, (j * 5 + 1) + 4), function(imp_j) {
          imputed_data_list[[exp]]$amelia_obj[[i]]$imputations[[imp_j]]
        }))
      })
      tmp_param_est <- sapply(tmp_amelia_obj_list, function(am_obj) {
        wa <- with.amelia(
          am_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,2]
      })
      tmp_sd_est <- sapply(tmp_amelia_obj_list, function(am_obj) {
        wa <- with.amelia(
          am_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,3]
      })
      tmp_param_est <- rowMeans(tmp_param_est)
      tmp_sd_est <- rowMeans(tmp_sd_est)
      model_list[[exp]][[i]] <- with.amelia(
        imputed_data_list[[exp]]$amelia_obj[[i]],
        lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
      )
      pooled_results[[exp]][[i]] <- summary(pool(as.mira(model_list[[exp]][[i]])))
      
      matrix <- matrix(
        c(tmp_param_est
          - qt(p = 0.975, df = df) * tmp_sd_est,
          tmp_param_est
          + qt(p = 0.975, df = df) * tmp_sd_est),
        ncol = 2
      )
      pooled_results[[exp]][[i]]$estimate <- tmp_param_est
      pooled_results[[exp]][[i]]$std.error <- tmp_sd_est
      
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


retransform_1col <- function(df, bn_obj, coln="x2") {
  df[[coln]] <- bestNormalize:::predict.bestNormalize(
    bn_obj[[coln]], newdata = df[[coln]], inverse = TRUE
  )
  df
}

retransform_columns <- function(imputed_data_list, columns = c("x2")) {
  for (exp in names(imputed_data_list)[-1]) {
    for (i in seq_along(imputed_data_list[[exp]]$data)) {
      for (col_index in seq_along(columns)) {
        imputed_data_list[[exp]]$data[[i]] <- retransform_1col(
          imputed_data_list[[exp]]$data[[i]],
          imputed_data_list[[exp]]$bn_objects[[i]][columns[col_index]],
          coln = columns[col_index]
        )
        for (imp_index in seq_along(imputed_data_list[[exp]]$amelia_obj[[i]]$imputations)) {
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



#########################################################################
#knn





get_graph_neighbors <- function(graph, var_name, labels) {
  edge_matrix <- str_split(
    names(graph@graph@edgeData@data), "\\|", simplify = TRUE
  )
  edge_matrix <- plyr::mapvalues(
    edge_matrix,
    from = as.character(seq_along(labels)),
    to = labels,
    warn_missing = FALSE
  )
  colnames(edge_matrix) <- c("from", "to")
  edges <- as_tibble(edge_matrix)
  edges$from <- str_replace(edges$from, "missing_", "")
  edges$to <- str_replace(edges$to, "missing_", "")
  edges <- edges %>%
    filter((from == var_name) | (to == var_name)) %>%
    as.matrix() %>%
    as.character() %>%
    unique()
  
  edges[edges != var_name]
}

impute_one_column_knn <- function(
    df, 
    df_preimputed, 
    col_name,
    feature_vector,
    nominal = FALSE,
    ...
) {
  extra_args = list(...)
  if ("k" %in% names(extra_args)) {
    k = extra_args[["k"]]
  } else {
    k = 2
  }
  if ("kernel" %in% names(extra_args)) {
    kernel = extra_args[["kernel"]]
  } else {
    kernel = "epanechnikov"
  }
  # cat("\nInside column name:", col_name)
  # cat("\nNominal:", nominal)
  # cat("\nFeature vec:", feature_vector)
  df_preimputed[[col_name]] <- df[[col_name]]
  missingness_index <- is.na(df[[col_name]])
  train_data <- df_preimputed[!missingness_index, ] %>%
    select(all_of(c(feature_vector, col_name)))
  test_data <- df_preimputed[missingness_index, ] %>%
    select(all_of(c(feature_vector, col_name)))
  if (nominal) {
    col_name_formula <- paste0("factor(", col_name, ")")
  } else {
    col_name_formula <- col_name
  }
  knn_imputations <- kknn(
    as.formula(paste0(col_name_formula, "~ .")),
    train_data, test_data, distance = 2,
    k = k,
    kernel = kernel
  )
  knn_imputations <- knn_imputations$fitted.values
  df[missingness_index, col_name] <- as.numeric(as.character(knn_imputations))
  df
}

impute_df_knn <- function(
    df, 
    df_preimputed, 
    mnar_cols,
    feature_vectors,
    nominal_columns = c("x6", "x7", "x8", "x9"),
    ...
) {
  for (col in colnames(df)) {
    if (col %in% mnar_cols) {
      nominal <- ifelse(col %in% nominal_columns, TRUE, FALSE)
      df <- impute_one_column_knn(
        df = df,
        df_preimputed = df_preimputed,
        feature_vector = feature_vectors[[col]],
        col_name = col,
        nominal = nominal,
        ...
      )
    } else {
      df[[col]] <- df_preimputed[[col]]  
    }
  }
  df
}

overall_impute_knn <- function(
    data_list, mnar_selection = "truth",
    graph_neighbors = TRUE, feature_vec_truth = NULL, 
    ...
) {
  for (exp in names(data_list)[-1]) {
    for (i in seq_along(data_list[[exp]]$data)) {
      # if a classified_data_list was given the detected mnar colmns are used
      if (is.list(mnar_selection)) {
        mnar_cols <- names(mnar_selection[[exp]]$detected_missingness_type[[i]][
          mnar_selection[[exp]]$detected_missingness_type[[i]] == "mnar"
        ])
      # for truth the ground truth missingness type is used
      } else if (mnar_selection == "truth") {
        mnar_cols <- data_list[[exp]]$missing[
          data_list[[exp]]$missing_type == "mnar"
        ]
      # in any other case the pre specified and fixed column selection is used
      } else {
        mnar_cols <- mnar_selection
      }
      feature_vectors <- sapply(mnar_cols, function(col) {
        if (is.list(mnar_selection) & graph_neighbors) {
          return(
            get_graph_neighbors(
              graph = mnar_selection[[exp]]$graph[[i]],
              var_name = col,
              labels = colnames(mnar_selection[[exp]]$data[[i]])
            )
          )
        } else {
          return(
            colnames(data_list[[exp]]$data[[i]])[
              !str_starts(colnames(data_list[[exp]]$data[[i]]), "missing_")
            ]
          )
        }
      }, simplify = FALSE, USE.NAMES = TRUE)
      
      if (!is.list(mnar_selection) && 
          ((mnar_selection == "truth") & graph_neighbors)) {
        feature_vectors <- feature_vec_truth
      }

      cat("\n\n-----------------------------")
      cat("\nExperiment:", exp, "Relative Frequency:", i)
      for (col in mnar_cols) {
        cat("\nMNAR Column:", col)
        cat("\nNeighbors:", str_flatten_comma(feature_vectors[[col]]))
      }
      
      data_list[[exp]]$knn_obj[[i]] <- impute_df_knn(
        df = data_list[[exp]]$data[[i]],
        df_preimputed = data_list[[exp]]$amelia_obj[[i]]$imputations$imp1,
        mnar_cols = mnar_cols,
        feature_vectors = feature_vectors,
        ...
      )
      data_list[[exp]]$knn_obj[[i]] <- sapply(
        data_list[[exp]]$amelia_obj[[i]]$imputations,
        FUN = function(imp) {
          impute_df_knn(
            df = data_list[[exp]]$data[[i]],
            df_preimputed = imp,
            mnar_cols = mnar_cols,
            feature_vectors = feature_vectors,
            ...
          )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
    }
  }
  data_list
}



evaluate_parameters_knn <- function(
    imputed_data_list,
    classified_data_list,
    form = as.formula("y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"),
    df = 990,
    estimated_missingness = TRUE
) {
  model_list <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  pooled_results <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  parameter_result <- list(mcar = vector("list", 3), mar = vector("list", 3), mnar = vector("list", 3))
  
  # raw data --------------------------------------------
  model_raw <- lm(
    formula = form,
    data = factorize_columns(imputed_data_list$raw$data[[1]])
  )
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
      imputation_length <- length(
        imputed_data_list[[exp]]$amelia_obj[[i]]$imputations
      )
      tmp_amelia_obj_list <- lapply((1:(imputation_length / 5)) - 1, function(j) {
        list(imputations = lapply(seq(j * 5 + 1, (j * 5 + 1) + 4), function(imp_j) {
          imputed_data_list[[exp]]$amelia_obj[[i]]$imputations[[imp_j]]
        }))
      })
      tmp_param_est <- sapply(tmp_amelia_obj_list, function(am_obj) {
        wa <- with.amelia(
          am_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,2]
      })
      tmp_sd_est <- sapply(tmp_amelia_obj_list, function(am_obj) {
        wa <- with.amelia(
          am_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,3]
      })
      tmp_param_est <- rowMeans(tmp_param_est)
      tmp_sd_est <- rowMeans(tmp_sd_est)
      model_list[[exp]][[i]] <- with.amelia(
        imputed_data_list[[exp]]$amelia_obj[[i]],
        lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
      )
      pooled_results[[exp]][[i]] <- summary(pool(as.mira(model_list[[exp]][[i]])))
      
      matrix <- matrix(
        c(tmp_param_est
          - qt(p = 0.975, df = df) * tmp_sd_est,
          tmp_param_est
          + qt(p = 0.975, df = df) * tmp_sd_est),
        ncol = 2
      )
      pooled_results[[exp]][[i]]$estimate <- tmp_param_est
      pooled_results[[exp]][[i]]$std.error <- tmp_sd_est
      
      pooled_results[[exp]][[i]] <- bind_cols(pooled_results[[exp]][[i]], matrix)
      colnames(pooled_results[[exp]][[i]])[7] <- "lower_bound"
      colnames(pooled_results[[exp]][[i]])[8] <- "upper_bound"
      
      # knn imputation ------------------------------------------------
      
      # if (estimated_missingness) {
      #   non_mnar_cols <- names(
      #     classified_data_list[[exp]]$detected_missingness_type[[i]][
      #       classified_data_list[[exp]]$detected_missingness_type[[i]] != "mnar"
      #     ]
      #   )
      # } else {
      #   non_mnar_cols <- names(
      #     classified_data_list[[exp]]$missing_type[
      #       classified_data_list[[exp]]$missing_type != "mnar"
      #     ]
      #   )
      # }
      
      tmp_knn_obj_list <- lapply((1:(imputation_length / 5)) - 1, function(j) {
        list(imputations = lapply(seq(j * 5 + 1, (j * 5 + 1) + 4), function(imp_j) {
          imputed_data_list[[exp]]$knn_obj[[i]][[imp_j]]
        }))
      })
      tmp_param_est <- sapply(tmp_knn_obj_list, function(knn_obj) {
        wa <- with.amelia(
          knn_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,2]
      })
      tmp_sd_est <- sapply(tmp_knn_obj_list, function(knn_obj) {
        wa <- with.amelia(
          knn_obj,
          lm(y ~ x1 + x2 + x3 + x4 + x5 + factor(x6) + x7 + x8 + x9) # CHANGE THIS MANUALLY
        )
        pr <- summary(pool(as.mira(wa)))
        pr[,3]
      })
      tmp_param_est <- rowMeans(tmp_param_est)
      tmp_sd_est <- rowMeans(tmp_sd_est)
      
      matrix_knn <- matrix(
        c(tmp_param_est
          - qt(p = 0.975, df = df) * tmp_sd_est,
          tmp_param_est
          + qt(p = 0.975, df = df) * tmp_sd_est),
        ncol = 2
      )
      
      knn_results <- tibble::tibble(
        est = tmp_param_est,
        coefficient = pooled_results[[exp]][[i]]$term,
        lower_bound = matrix_knn[, 1],
        upper_bound = matrix_knn[, 2]
      ) %>%
        mutate(type = "knn_imputed")
      
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
        bind_rows(knn_results)
      
      parameter_result[[exp]][[i]] <- parameter_result[[exp]][[i]] %>%
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