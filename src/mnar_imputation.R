library(kknn)
load("data/visualized_mixed_data.RData")

impute_one_column_knn <- function(
    df, 
    df_preimputed, 
    col_name,
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
  df_preimputed[[col_name]] <- df[[col_name]]
  missingness_index <- is.na(df[[col_name]])
  train_data <- df_preimputed[!missingness_index, ]
  test_data <- df_preimputed[missingness_index, ]
  knn_imputations <- kknn(
    as.formula(paste0(col_name, "~.")),
    train_data, test_data, distance = 2,
    k = k,
    kernel = kernel
  )$fitted.values
  df[missingness_index, col_name] <- knn_imputations
  df
}

impute_df_knn <- function(df, df_preimputed, mnar_cols, ...) {
  for (col in colnames(df)) {
    if (col %in% mnar_cols) {
      df <- impute_one_column_knn(
        df = df,
        df_preimputed = df_preimputed,
        col_name = col,
        ...
      )
    } else {
      df[[col]] <- df_preimputed[[col]]  
    }
  }
  df
}

overall_impute_knn <- function(data_list) {
  for (i in seq_along(data_list$mnar$data)) {
    data_list$mnar$knn_obj[[i]] <- impute_df_knn(data_list$mnar$data[[i]],
                                                   data_list$mnar$amelia_obj[[i]]$imputations$imp1,
                                                   data_list$mnar$missing[data_list$mnar$missing_type == "mnar"])
  }
  data_list
}

knn_imputed_mixed_data_list <- overall_impute_knn(visualized_mixed_data_list)

# selbe diagnostische plots
knn_density_comparison_plots <- function(knn_imputed_data_list, col_name, rel_miss) {
  rel_miss_ind <- which(
    rel_miss == knn_imputed_data_list$mnar$relative_missingness
  )
  if (length(rel_miss_ind) == 0) {
    stop("This Relative missingness does not exist.")
  }
  knn_imputed_data_list$mnar$dens_plots[[col_name]][[rel_miss_ind]] +
    geom_line(
      data = data.frame(
        x = density(
          knn_imputed_data_list$mnar$knn_obj[[rel_miss_ind]][[col_name]]
        )$x,
        y = density(
          knn_imputed_data_list$mnar$knn_obj[[rel_miss_ind]][[col_name]]
        )$y
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

#scheint für kleine k besser zu funktionieren
# x4 als input funktioniert nicht weil x4 nicht missing ist!
# funktioniert nicht für categorical variablen aber bei categorical lieber 
# kreuztabelle machen weil man bei Balken e nichts sieht
# x1 und x7 sind mcar. imputation dazu ist amelia sieht aber anders aus als rote linie
# weil wir erste imputation nehmen und nicht mitteln
knn_density_comparison_plots(
  knn_imputed_data_list = knn_imputed_mixed_data_list,
  col_name = "x2",
  rel_miss = 0.3
)

