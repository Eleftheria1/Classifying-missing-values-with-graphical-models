# load requiered packages
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


# overall_impute: Performs multiple imputation using Amelia II for missing data in the datasets of data_list.
# Parameters:
#   - data_list: A list containing multiple datasets for different experiments MCAR, MAR, MNAR.
#   - nominals: A vector of column names indicating nominal/categorical variables (default = NULL).
#   - m: The number of imputed datasets to generate (default = 5).
#   - p2s: Indicates the screen output.
#   - sqrts: A vector of column names for which the square root transformation is applied before imputation (default = NULL).
# Description:
# The 'overall_impute' function applies multiple imputation using the Amelia II package to impute missing data
# in the datasets of 'data_list'. For each dataset in each experiment, it uses the 'amelia' function to perform
# imputation.
# The function updates 'data_list' by storing the imputation results in the 'amelia_obj' field for each dataset.
# Returns the updated data_list containing the Amelia II imputation objects.
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

# factorize_columns: Converts specified columns in a data frame to factors if they have more than one unique value.
# Parameters:
#   - df: Data frame containing the columns to be factorized.
#   - columns: A character vector specifying the column names to be factorized (default = "x6").
# Description:
# For each column specified in 'columns', it checks the number of unique values using 'length(unique(df[[col]]))'.
# If the number of unique values is greater than one, it converts the column to a factor using the 'factor' function.
# The function then returns the updated data frame 'df' with the factorized columns.
factorize_columns <- function(df, columns = c("x6")) {
  for (col in columns) {
    if (length(unique(df[[col]])) > 1) {
      df[[col]] <- factor(df[[col]])
    }
  }
  df
}

# evaluate_parameters: Evaluates the parameters (coefficients) of the linear model for different missingness types (MCAR, MAR, MNAR)
# using imputed datasets and complete case data.
# Parameters:
#   - imputed_data_list: A list containing imputed datasets for different missingness types.
#   - form: The formula representing the linear model to be evaluated (default formula: y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9).
#   - df: Degrees of freedom for confidence intervals (default = 990).
# Description:
# The 'evaluate_parameters' function performs parameter evaluation for different missingness types (MCAR, MAR, MNAR) using imputed datasets
# and complete case data. It initializes three lists 'model_list', 'pooled_results', and 'parameter_result' to store the results
# for each missingness type.
# The function first fits the linear model ('model_raw') to the raw data using the specified formula ('form').
# It computes the model's coefficients and confidence intervals at a 95% level of confidence and stores the results in 'raw_results' tibble.
# The function then iterates through each experiment and imputed dataset to fit the model.
# It calculates the pooled estimates and standard errors for the model parameters using the 'summary' and 'pool' functions.
# The pooled results are stored in 'pooled_results' list. Additionally, the function computes complete case results by fitting
# the linear model to the data with missing values removed.
# The final 'parameter_result' list combines the results for imputed data, complete case data, and raw data, along with the type of each result.
# Returns a list containing the 'model_list', 'pooled_results', and 'parameter_result' with the parameter evaluation results
# for different missingness types.
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

# normalize_1col: Normalizes a single column in a data frame using the 'bestNormalize' function from the 'bestNormalize' package.
# Parameters:
#   - df: Data frame containing the column to be normalized.
#   - col: Name of the column to be normalized (default = "x2").
# Description:
# The 'normalize_1col' function calculates the best transformation for normalization and applies it to the column. 
# The 'standardize' parameter is set to FALSE to preserve the original scale
# of the data. The normalized column replaces the original column in the data frame.
# The function returns a list containing the 'bn_obj', which contains information about the transformation applied by 'bestNormalize',
# and the updated data frame with the normalized column.
normalize_1col <- function(df, col = "x2") {
  bn_obj <- bestNormalize(df[[col]], warn = TRUE, standardize = FALSE)
  df[[col]] <- bn_obj$x.t
  list(bn_obj, df)
}

# normalize_columns: Normalizes specified columns in datasets within the data_list using 'bestNormalize' function.
# Parameters:
#   - data_list: A list containing multiple datasets for different experiments.
#   - columns: A character vector specifying the column names to be normalized (default = "x2").
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

# retransform_1col: Retransforms a single column in a data frame using the 'predict.bestNormalize' function from the 'bestNormalize' package.
# Parameters:
#   - df: Data frame containing the column to be retransformed.
#   - bn_obj: List containing the transformation information obtained from 'bestNormalize' for the specified column.
#   - coln: Name of the column to be retransformed (default = "x2").
# Description:
# The 'retransform_1col' function uses the transformation information stored in 'bn_obj' to perform the inverse transformation on the column.
# The 'inverse' parameter is set to TRUE to indicate that the retransformation should be performed in the reverse direction.
# The retransformed column replaces the original column in the data frame.
# The function returns the updated data frame with the retransformed column.
retransform_1col <- function(df, bn_obj, coln="x2") {
  df[[coln]] <- bestNormalize:::predict.bestNormalize(
    bn_obj[[coln]], newdata = df[[coln]], inverse = TRUE
  )
  df
}

# retransform_columns: Retransforms specified columns in datasets within the imputed_data_list using 'predict.bestNormalize' function.
# Parameters:
#   - imputed_data_list: A list containing multiple imputed datasets for different experiments MCAR, MAR, MNAR.
#   - columns: A character vector specifying the column names to be retransformed (default = "x2").
# Description:
# The 'retransform_columns' function iterates through each experiment and imputed dataset to retransform the columns individually.
# For each column in the specified 'columns', it calls the 'retransform_1col' function to apply retransformation using 'predict.bestNormalize'.
# The retransformed columns replace the normalized columns in the datasets.
# The function returns the updated imputed_data_list with the retransformed columns.
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




# get_graph_neighbors: Get the neighboring variables of a specified variable in a given graph.
# Parameters:
#   - graph: An igraph object representing the graph structure.
#   - var_name: The name of the variable for which neighbors need to be found.
#   - labels: A character vector containing the variable names (node labels) in the graph.
# Description:
# The 'get_graph_neighbors' function extracts the neighbors of a specified variable ('var_name') in the graph.
# It does this by processing the edge matrix of the graph.
# The function converts the edge matrix into a data frame ('edges') with columns "from" and "to" representing the connected nodes.
# Next, the function filters the data frame to include only the rows where the specified variable ('var_name') appears either as the "from" or "to" variable.
# The resulting data frame is then converted to a character matrix and unique values are extracted to obtain the neighboring variable names.
# The function returns the unique neighboring variable names excluding the specified variable ('var_name').
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


# impute_one_column_knn: Impute missing values in a single column using K-nearest neighbor imputation.
# Parameters:
#   - df: Data frame containing the column to be imputed.
#   - df_preimputed: Data frame with the amelia imputations.
#   - col_name: Name of the column to be imputed.
#   - feature_vector: A character vector containing the names of the feature (predictor) variables to be used in the imputation.
#   - nominal: Logical value indicating if the column contains nominal/categorical data (default = FALSE).
#   - ...: Additional arguments that can be passed to the 'kknn' function.
# Description:
# The 'impute_one_column_knn' function performs K-nearest neighbor imputation for the specified column ('col_name') in the data frame ('df').
# The function identifies the missing values in the column and creates training and test data sets accordingly.
# The imputation is carried out using the 'kknn' function from the 'kknn' package.
# The imputed values are replaced in the original data frame ('df') at the positions of the missing values.
# The function returns the updated data frame with the imputed column.
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

# impute_df_knn: Impute missing values in multiple columns of a data frame using K-nearest neighbor imputation.
# Parameters:
#   - df: Data frame with missing values to be imputed.
#   - df_preimputed: Data frame containing the Amelia imputations.
#   - mnar_cols: A character vector specifying the column names with missing not at random (MNAR) data.
#   - feature_vectors: A list of feature vectors containing the names of predictor variables for imputation.
#   - nominal_columns: A character vector specifying the column names containing nominal/categorical data (default = c("x6", "x7", "x8", "x9")).
#   - ...: Additional arguments that can be passed to the 'kknn' function for each imputation.
# Description:
# The 'impute_df_knn' function applies K-nearest neighbor imputation to multiple columns in the data frame ('df').
# It iterates through each column in 'df', and for each MNAR column (specified in 'mnar_cols'), it calls 'impute_one_column_knn'.
# The 'feature_vectors' list should have a feature vector (list of predictor variable names) for each column to be imputed.
# If a column is not in 'mnar_cols', the Amelia imputations from 'df_preimputed' are kept in the imputed data frame.
# The 'nominal_columns' parameter allows the user to specify the columns that contain nominal/categorical data for proper imputation.
# The function returns the data frame with the missing values imputed using K-nearest neighbor imputation.
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

# overall_impute_knn: Perform overall imputation using K-nearest neighbor imputation for MNAR (missing not at random) data.
# Parameters:
#   - data_list: A list of data frames representing different experiments MCAR, MAR, MNAR.
#   - mnar_selection: A parameter specifying how to select MNAR columns for imputation. It can be one of the following:
#                     - "truth" (default): Use the true MNAR columns..
#                     - A character vector of column names: Use the provided fixed selection of MNAR columns.
#                     - A list of classification results: Use the detected MNAR columns from the classification results.
#   - graph_neighbors: A logical value indicating whether to use graph neighbors for feature vectors (default = TRUE).
#   - feature_vec_truth: A list of feature vectors corresponding to the ground truth missingness type (used when mnar_selection is "truth" and graph_neighbors is TRUE).
#   - ...: Additional arguments to be passed to 'impute_df_knn' function for imputation.
# Description:
# The 'overall_impute_knn' function iterates through each experiment and relative missigness in the 'data_list'.
# The function determines the MNAR columns to be imputed based on 'mnar_selection'.
# If 'graph_neighbors' is TRUE and 'mnar_selection' is a list of classification results, it uses graph neighbors to form feature vectors for KNN imputation.
# The function then applies K-nearest neighbor imputation using 'impute_df_knn' to each MNAR column in each data frame of the experiment.
# The results of imputation are stored in 'data_list$knn_obj'.
# The function returns the updated 'data_list' with the imputed values using K-nearest neighbor imputation.
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


# evaluate_parameters_knn: Evaluate parameters using K-nearest neighbor (KNN) imputation for MNAR (missing not at random) data.
# Parameters:
#   - imputed_data_list: A list of data frames representing different experiments MCAR, MAR, MNAR, each containing imputed data.
#   - classified_data_list: A list of data frames representing different experiments MCAR, MAR, MNAR, each containing classified missingness types.
#   - form: A formula specifying the linear regression model to be used for parameter evaluation (default formula: "y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9").
#   - df: Degrees of freedom for confidence intervals (default = 990).
#   - estimated_missingness: A logical value indicating whether to use estimated missingness type (default = TRUE).
# Description:
# The 'evaluate_parameters_knn' function evaluates parameters using K-nearest neighbor (KNN) imputation for MNAR data.
# It calculates the parameters for each experiment and relative missigness in the 'imputed_data_list'.
# The MNAR columns to be imputed are determined based on the 'classified_data_list'.
# The function then applies K-nearest neighbor imputation using 'impute_df_knn' to each MNAR column in each data frame of the experiment.
# The KNN-imputed data is used to fit linear regression models, and parameter estimates are obtained for each imputation.
# The function also calculates confidence intervals for the parameter estimates based on 'df'.
# The results of parameter estimates and confidence intervals are stored in 'parameter_result' list.
# The function returns 'model_list', 'pooled_results', and 'parameter_results' containing the respective model objects,
# pooled results, and parameter estimates with confidence intervals.
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