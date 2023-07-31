# Get required packages
library(tidyverse)

use_mvpc <- TRUE
if (use_mvpc) {
  source("mvpc/CITest.R")
  source("mvpc/MissingValuePC.R")
}


# Visualization Utils

# plot_graph_fit: Plots a graph representing a fitted m-graph using igraph.
# Parameters:
#   - graph_fit: A fitted m-graph object.
#   - labels: A vector of labels for the graph vertices.
#   - experiment: Experiment identifier for the graph plot e.g. MCAR/MAR or MNAR.
#   - rel_missingess: Relative missingness value for the graph plot e.g. 0.1, 0.3 or 0.6.
# The function converts the 'graph_fit' object to an igraph object and sets vertex labels.
# Then, it plots the igraph object with specified aesthetics
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

# fit_mvpc_graph: Fits a mvpc-graph using given data and the code of Ruibo Tu 
# (Tu, Ruibo ; Zhang, Cheng ; Ackermann, Paul ; Mohan, Karthika ; Kjellström,
# Hedvig ; Zhang, Kun: Causal discovery in the presence of missing data. In: The
# 22nd International Conference on Artificial Intelligence and Statistics PMLR, 2019)
# Parameters:
#   - data: Data matrix with variables in columns and samples in rows.
#   - relative_missingness: Relative missingness value for the experiment e.g. 0.1, 0.3 or 0.6.
#   - experiment: Experiment identifier e.g. MCAR/MAR or MNAR.
#   - alpha: Significance level for partial correlation tests (default = 0.01).
#   - indep_test: Function for independence testing (default = gaussCItest.td).
#   - corr_method: Function for computing partial correlations (default = gaussCItest.permc).
#   - plot_graph: Boolean; if TRUE, plots the fitted mvpc graph using plot_graph_fit function.
# The function fits the mvpc graph using the 'mvpc' function with specified parameters:
#   - 'suffStat' as a list with the data,
#   - 'indepTest' for independence testing function,
#   - 'corrMethod' for partial correlation computation function,
#   - 'p' as the number of columns in 'data',
#   - and 'alpha' as the significance level for the tests.
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

# fit_pc_graph: Fits a pc graph using given data and plots the graph if requested.
# Parameters:
#   - data: Data matrix with variables in columns and samples in rows.
#   - relative_missingness: Relative missingness (e.g. 0.1, 0.3, 0.6) value for the experiment.
#   - experiment: Experiment identifier (e.g. MCAR, MAR, MNAR).
#   - alpha: Significance level for independence tests (default = 0.01).
#   - plot_graph: Boolean; if TRUE, plots the fitted PC graph using plot_graph_fit function.
# Returns the fitted PC graph object.
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

# apply_graph_fitting: Fits a graph of specified type (pc or mvpc) to a dataset in the data_list.
# Parameters:
#   - data_list: A list containing multiple datasets for different experiments MCAR,MAR,MNAR.
#   - experiment: Experiment identifier for the dataset to be analyzed.
#   - rel_missingness: Relative missingness value for the dataset to be analyzed (e.g. 0.1, 0.3, 0.6).
#   - type: Type of graph to fit ("pc" or "mvpc").
#   - alpha: Significance level for independence tests (default = 0.01).
#   - plot_graph: Boolean; if TRUE, plots the fitted graph.
# The function first matches the 'type' argument to "pc" or "mvpc". It then finds the index of the dataset
# corresponding to the specified 'experiment' and 'rel_missingness' in 'data_list'. 
# Depending on the 'type', it calls 'fit_pc_graph' or 'fit_mvpc_graph' to fit the corresponding graph
# to the dataset and optionally plots the fitted graph using 'plot_graph_fit' function.
# Returns the fitted graph object.
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


# fit_and_visualize_graphs: Fits and visualizes graphs (pc or mvpc) to datasets in data_list.
# Parameters:
#   - data_list: A list containing multiple datasets for different experiments (MCAR, MAR, MNAR).
#   - type: Type of graph to fit ("pc" or "mvpc").
#   - alpha: Significance level for independence tests (default = 0.01).
#   - plot_graph: Boolean; if TRUE, plots the fitted graphs.
# Description:
# The 'fit_and_visualize_graphs' function iterates over datasets in 'data_list' for each experiment
# and relative missingness level. Depending on the specified 'type', it calls 'fit_pc_graph' or 'fit_mvpc_graph'
# to fit the corresponding graph to the dataset. If 'plot_graph' is TRUE, it plots the fitted graph using
# 'plot_graph_fit' function.
# The function updates 'data_list' by adding the fitted graph objects to the 'graph' field for each dataset.
# Returns the updated data_list containing the fitted graph objects.
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

# classify_missingness_variable: Classifies the type of missingness (MCAR, MAR, or MNAR) for a given variable
# based on the graph structure and its relationships with other variables.
# Parameters:
#   - graph: The fitted graph object, representing the conditional dependencies between variables.
#   - missingness_variable: The variable for which the missingness type needs to be determined.
#   - all_missing_variables: A vector containing all variables with potential missingness.
#   - labels: A vector of labels representing variable names.
# Description:
# The 'classify_missingness_variable' function analyzes the 'graph' structure to classify the missingness
# type for the specified 'missingness_variable'. It extracts the edges from the 'graph' and compares them with
# the provided 'labels'. The edges are then categorized as either directed or undirected relationships.
# The function identifies the presence of missingness variables in the graph edges to determine if the missingness
# is MCAR, MAR, or MNAR.
# Returns a character indicating the missingness type: "mcar" (missing completely at random),
# "mar" (missing at random), or "mnar" (missing not at random).
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
  if (purrr::is_empty(edge_matrix)) return("mcar")
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

# classify_all_missingness_variables: Classifies the missingness type (MCAR, MAR, or MNAR) for all specified variables
# based on the graph structure and their relationships with other variables.
# Parameters:
#   - graph: The fitted graph object, representing the conditional dependencies between variables.
#   - all_missingness_variables: A vector containing all variables with potential missingness.
#   - labels: A vector of labels representing variable names.
# Description:
# The 'classify_all_missingness_variables' function applies 'classify_missingness_variable' to each variable
# in 'all_missingness_variables' to determine the missingness type for each variable. For each variable, it analyzes
# the 'graph' structure to classify its missingness as MCAR, MAR, or MNAR.
# Returns a vector of characters indicating the missingness type for each specified variable.
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

# add_classification: Classifies missingness type (MCAR, MAR, or MNAR) for all datasets in data_list
# and adds the classification results to the data_list.
# Parameters:
#   - data_list: A list containing multiple datasets for different experiments (MCAR, MAR, or MNAR).
# Description:
# The 'add_classification' function iterates over datasets in 'data_list' for each experiment
# and relative missingness level. For each dataset, it calls 'classify_all_missingness_variables'
# to determine the missingness type for all specified missingness variables based on the fitted graph.
# The missingness types are then added to the 'detected_missingness_type' field of each dataset
# in the data_list.
# Returns the updated data_list containing the added missingness classification information.
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
