#############################################################
#                  Fit the graphical models
#############################################################

# Get required packages & data
library(tidyverse)
load("data/multi_normal.RData")
load("data/mixed_data.RData")
source("mvpc/CITest.R")
source("mvpc/MissingValuePC.R")


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

apply_mvpc <- function(
  data_list,
  experiment,
  rel_missingness,
  alpha = 0.01,
  plot_graph = TRUE
  ) {
  rel_missingness_index <- which(
    data_list[[experiment]]$relative_missingness == rel_missingness
  )
  fit_mvpc_graph(
    data = as.data.frame(data_list[[experiment]]$data[[rel_missingness_index]]),
    relative_missingness = rel_missingness,
    experiment = experiment,
    plot_graph = plot_graph,
    alpha = alpha
  )
}
set.seed(123)
fitted_mvpc_graph <- apply_mvpc(
  multi_normal_data_list,
  experiment = "mnar_x2",
  rel_missingness = 0.3,
  alpha = 0.01
)
fitted_mvpc_graph <- apply_mvpc(
  mixed_data_list,
  experiment = "mar",
  rel_missingness = 0.1,
  alpha = 0.01
)

# fit and visualize all at the same time
for (experiment in names(multi_normal_data_list)[2:length(multi_normal_data_list)]) {
  for (rel_missing_i in seq_along(multi_normal_data_list[[experiment]]$relative_missingness)) {
    multi_normal_data_list[[experiment]]$graph[[rel_missing_i]] <- fit_mvpc_graph(
      data = multi_normal_data_list[[experiment]]$data[[rel_missing_i]],
      relative_missingness = multi_normal_data_list[[experiment]]$relative_missingness[[rel_missing_i]],
      experiment = experiment,
      plot_graph = TRUE, alpha = 0.01
    )
  }
}

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

# classify_missingness_variable(
#   fitted_mvpc_graph,
#   missingness_variable = "x5",
#   all_missing_variables = c("x2", "x5"),
#   labels = c("x1", "x2", "x3", "x4", "x5", "y", "missing_x2", "missing_x5")
# )

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

# detect missingness for the full list
for (experiment in names(multi_normal_data_list)[2:length(multi_normal_data_list)]) {
  for (rel_missing_i in seq_along(multi_normal_data_list[[experiment]]$relative_missingness)) {
    multi_normal_data_list[[experiment]]$detected_missingness_type[[rel_missing_i]] <- classify_all_missingness_variables(
      graph = multi_normal_data_list[[experiment]]$graph[[rel_missing_i]],
      all_missingness_variables = multi_normal_data_list[[experiment]]$missing, 
      labels = colnames(multi_normal_data_list[[experiment]]$data[[rel_missing_i]])
    )
  }
}


classification_results <- data.frame(
  experiment = rep("", 3 * 3),
  relative_missingness = numeric(3 * 3),
  correct_classification = rep(FALSE, 3 * 3)
)
index <- 1
for (experiment in names(multi_normal_data_list)[2:length(multi_normal_data_list)]) {
  for (rel_missing_i in seq_along(multi_normal_data_list[[experiment]]$relative_missingness)) {
    classification_results$experiment[index] <- experiment
    classification_results$relative_missingness[index] <- multi_normal_data_list[[experiment]]$relative_missingness[rel_missing_i]
    classification_results$correct_classification[index] <- all(
      sapply(multi_normal_data_list[[experiment]]$detected_missingness_type, function(detection) {
        print(detection)
        print(multi_normal_data_list[[experiment]]$missing_type)
        all(detection == multi_normal_data_list[[experiment]]$missing_type)
      }, simplify = TRUE)
    )
    index <- index + 1
  }
}

classification_results %>%
  mutate(
    relative_missingness = as.factor(relative_missingness),
    correct_classification = as.character(correct_classification),
    experiment = str_to_upper(str_replace(experiment, "_", " "))
  ) %>%
  ggplot(
    aes(x = relative_missingness, y = experiment, fill = correct_classification)
  ) +
  geom_tile(col = "white") +
  scale_fill_manual(
    values = c("FALSE" = "orange", "TRUE" = "darkgreen"),
    name = "Classification",
    labels = c("TRUE" = "Correct", "FALSE" = "Wrong")
  ) +
  labs(x = "Relative Missingness", y = "Experiment") +
  theme_minimal()






#####################################################################
# OLD
#####################################################################






#BiocManager::install("graph")
#BiocManager::install("RBGL")
#BiocManager::install("Rgraphviz")
library(graph)
library(igraph)
library(gRbase)
library(Rgraphviz)
library(RBGL)

library(ggm)
library(pcalg)

?ggm
complete_data <- as.data.frame(multi_normal_data_list$raw[[1]])[,-1]

graph <- DAG(complete_data$y ~ complete_data$x1 + complete_data$x2 + 
                      complete_data$x3 + complete_data$x4)
plot(as(graph, "graphNEL"))
graph_object <- as(graph, "graphNEL")
nodes(graph_object) <- as.character(nodes(graph_object))
edges <- buildEdgeList(graph_object)
plot(graph_object)
#fitted_graph <- fitDag(graph_object, S.carc)






suffStat <- list(C = cor(complete_data), n = nrow(complete_data))
pc_algorithm <- pc(suffStat, indepTest = gaussCItest,
             p = ncol(complete_data), alpha = 0.01)
stopifnot(require(Rgraphviz))
plot(graph_object, main = "") ; plot(pc_algorithm, main = "")


missing_data <- as.data.frame(multi_normal_data_list$mar_x2[[1]][[2]])
correlation <- cor(missing_data)
cor_complete_cases_only <- cor(na.omit(missing_data))
correlation[2,] <- cor_complete_cases_only[2,]
correlation[,2] <- cor_complete_cases_only[,2]
suffStat_missing <- list(C = unname(correlation), n = nrow(missing_data[, 1:5]), data = missing_data[, 1:5])
pc_algorithm_missing <- pc(suffStat_missing, indepTest = gaussCItest,
                   p = ncol(missing_data), alpha = 0.01)
plot(pc_algorithm_missing, main = "")
suff_missing <- list(data = missing_data)

mvpc_mar <- mvpc(suffStat = suff_missing, indepTest =  gaussCItest.td,
                 corrMethod =  gaussCItest.permc,
                 p = 7, alpha = 0.01)
plot(mvpc_mar)

mnar_data <- as.data.frame(multi_normal_data_list$mnar_x2$data[[1]])
suff_missing_mnar <- list(data = mnar_data)
mvpc_mnar <- mvpc(suffStat = suff_missing_mnar, indepTest =  gaussCItest.td,
                 corrMethod =  gaussCItest.permc,
                 p = 7, alpha = 0.01)
plot(mvpc_mnar)


correlation_mnar <- cor(mnar_data)
cor_complete_cases_only_mnar <- cor(na.omit(mnar_data))
correlation_mnar[2,] <- cor_complete_cases_only_mnar[2,]
correlation_mnar[,2] <- cor_complete_cases_only_mnar[,2]
correlation_mnar[5,] <- cor_complete_cases_only_mnar[5,]
correlation_mnar[,5] <- cor_complete_cases_only_mnar[,5]
a <- mnar_data %>% select(-x2)
a_cor <- cor(na.omit(a))
correlation_mnar[5,7] <- a_cor[4,6]
correlation_mnar[7,5] <- a_cor[6,4]

set.seed(123)
suffStat_mnar <- list(C = correlation_mnar, n = nrow(mnar_data))
suffStat_mnar <- list(C = cor(mnar_data, use = "pairwise.complete.obs"), n = nrow(mnar_data))
pc_algorithm_mnar <- pc(suffStat_mnar, indepTest = gaussCItest,
                           p = ncol(mnar_data), alpha = 0.01)


plot(pc_algorithm_mnar, main = "")

mcar_data <- as.data.frame(multi_normal_data_list$mcar_x2[[1]][[1]])
suffStat_mcar <- list(C = cor(mcar_data), n = nrow(mcar_data))
pc_algorithm_mcar <- pc(suffStat_mcar, indepTest = gaussCItest,
                        p = ncol(mcar_data), alpha = 0.01)

plot(pc_algorithm_mcar, main = "")

library(gRim)

# 
# library(gRbase)
# data(milkcomp1, package='gRbase')
# SS <- CGstats(milkcomp1, varnames=c("treat","fat","protein",
#                                     "lactose"))
# SS
# can.parms<-CGstats2mmodParms(SS,type="ghk")
# 
# milkmod <- mmod(~treat*fat*protein + protein*lactose, data=milkcomp1)
# str(milkmod)
# 
# plot(milkmod)
# 
# glist <- ~treat:fat:protein+protein:lactose
# milk <- mmod(glist, data=milkcomp1)
# summary(milk)
# coef(milk, type="ghk")
# mm <- mmod(~.^., data=milkcomp1)
# mm2 <- stepwise(mm, k=log(nrow(milkcomp1)), details=0)
# plot(mm2)
# 
# complete_data$missing <- as.factor(complete_data$missing)
# mm <- mmod(~.^., data=complete_data)
# mm2 <- stepwise(mm, k=log(nrow(complete_data)), details=0)
# plot(mm2)
# 
# missing_data$missing <- as.factor(missing_data$missing)
# mm <- mmod(~.^., data=missing_data)
# mm2 <- stepwise(mm, k=log(nrow(missing_data)), details=0)
# plot(mm2)











