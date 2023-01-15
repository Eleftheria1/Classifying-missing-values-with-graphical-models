# Get required packages
library(tidyverse)

use_mvpc <- TRUE
if (use_mvpc) {
  source("mvpc/CITest.R")
  source("mvpc/MissingValuePC.R")
}


# Visualization Utils
plot_graph_fit <- function(graph_fit, labels, experiment, rel_missingess) {
  igraph_obj <- igraph::graph_from_graphnel(graph_fit@graph)
  igraph_obj <- igraph::set.vertex.attribute(
    igraph_obj,
    "label",
    value = labels
  )
  igraph::plot.igraph(
    igraph_obj,
    mark.col = "blue",
    vertex.color = "lightblue",
    vertex.size = 20,
    edge.arrow.size = 0.5,
    label.cex = 12,
    main = paste(
      "m-graph for experiment:", experiment,
      "and relative missingness:", rel_missingess
    )
  )
  # iplotPC(graph_fit, labels = labels)
}

# Graph Fitting

fit_mvpc_graph <- function(
    data, relative_missingness, experiment,
    alpha = 0.01,
    indep_test = gaussCItest.td, 
    corr_method = gaussCItest.permc,
    plot_graph = TRUE
) {
  mvpc_fit <- mvpc(
    suffStat = list(data = as.data.frame(data)),
    indepTest = indep_test,
    corrMethod =  corr_method,
    p = ncol(data),
    alpha = alpha
  )
  if (plot_graph) {
    plot_graph_fit(
      mvpc_fit,
      labels = colnames(data),
      experiment = experiment,
      rel_missingess = relative_missingness
    )
  }
  mvpc_fit
}

fit_pc_graph <- function(
    data, relative_missingness, experiment,
    alpha = 0.01,
    plot_graph = TRUE
) {
  filled_cor <- cor(as.data.frame(data), use = "pairwise.complete.obs")
  # set NAs to zero according to the assumption of no self masked missingness
  filled_cor[is.na(filled_cor)] <- 0
  pc_fit <- pcalg::pc(
    suffStat = list(
      C = filled_cor,
      n = nrow(data)
    ),
    indepTest = pcalg::gaussCItest,
    p = ncol(data),
    alpha = alpha
  )
  if (plot_graph) {
    plot_graph_fit(
      pc_fit,
      labels = colnames(data),
      experiment = experiment,
      rel_missingess = relative_missingness
    )
  }
  pc_fit
}

apply_graph_fitting <- function(
    data_list,
    experiment,
    rel_missingness,
    type = c("pc", "mvpc"),
    alpha = 0.01,
    plot_graph = TRUE
) {
  type = match.arg(type)
  rel_missingness_index <- which(
    data_list[[experiment]]$relative_missingness == rel_missingness
  )
  if (type == "mvpc") {
    fit_mvpc_graph(
      data = as.data.frame(data_list[[experiment]]$data[[rel_missingness_index]]),
      relative_missingness = rel_missingness,
      experiment = experiment,
      plot_graph = plot_graph,
      alpha = alpha
    )
  } else if (type == "pc") {
    fit_pc_graph(
      data = as.data.frame(data_list[[experiment]]$data[[rel_missingness_index]]),
      relative_missingness = rel_missingness,
      experiment = experiment,
      plot_graph = plot_graph,
      alpha = alpha
    )
  } else {
    stop("Type not implemented")
  }
}

fit_and_visualize_graphs <- function(
    data_list,
    type = c("pc", "mvpc"),
    alpha = 0.01,
    plot_graph = FALSE
) {
  type = match.arg(type)
  for (experiment in names(data_list)[2:length(data_list)]) {
    for (rel_missing_i in seq_along(data_list[[experiment]]$relative_missingness)) {
      if (type == "mvpc") {
        data_list[[experiment]]$graph[[rel_missing_i]] <- fit_mvpc_graph(
          data = data_list[[experiment]]$data[[rel_missing_i]],
          relative_missingness = data_list[[experiment]]$relative_missingness[[rel_missing_i]],
          experiment = experiment,
          plot_graph = plot_graph, alpha = alpha
        )
      } else {
        data_list[[experiment]]$graph[[rel_missing_i]] <- fit_pc_graph(
          data = data_list[[experiment]]$data[[rel_missing_i]],
          relative_missingness = data_list[[experiment]]$relative_missingness[[rel_missing_i]],
          experiment = experiment,
          plot_graph = plot_graph, alpha = alpha
        )
      }
    }
  }
  data_list
}

# Graph based missingness classification
classify_missingness_variable <- function(
    graph,
    missingness_variable,
    all_missing_variables,
    labels
) {
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
  
  missingness_ind <- paste0("missing_", missingness_variable)
  # detect MCAR
  # no edge to and from the missingness indicator allowed
  if (all(missingness_ind != edge_matrix)) {
    return("mcar")
  }
  
  # detect MAR
  edge_matrix <- edge_matrix[
    rowSums(missingness_ind != edge_matrix) != 2,
    ,
    drop = FALSE
  ]
  n_edges <- nrow(edge_matrix)
  for (row in 1:n_edges) {
    if (any(edge_matrix[row, ] %in% all_missing_variables)) {
      return("mnar")
    }
    if (nrow(edge_matrix) > 1) {
      inverse_edge <- c(edge_matrix[row, 2], edge_matrix[row, 1])
      match_inverse <- sapply(
        seq(n_edges)[seq(n_edges) != row],
        function(other_row) {
          all(edge_matrix[other_row, ] == inverse_edge)
        }
      )
      if (any(match_inverse)) {
        return("mnar")
      }
    }
  }
  "mar"
}

classify_all_missingness_variables <- function(
    graph,
    all_missingness_variables,
    labels
) {
  sapply(all_missingness_variables, function(missingness_variable) {
    classify_missingness_variable(
      graph = graph,
      missingness_variable = missingness_variable,
      all_missing_variables = all_missingness_variables,
      labels = labels
    )
  }, simplify = TRUE)
}

add_classification <- function(data_list) {
  for (experiment in names(data_list)[2:length(data_list)]) {
    for (rel_missing_i in seq_along(data_list[[experiment]]$relative_missingness)) {
      data_list[[experiment]]$detected_missingness_type[[rel_missing_i]] <- classify_all_missingness_variables(
        graph = data_list[[experiment]]$graph[[rel_missing_i]],
        all_missingness_variables = data_list[[experiment]]$missing, 
        labels = colnames(data_list[[experiment]]$data[[rel_missing_i]])
      )
    }
  }
  data_list
}
