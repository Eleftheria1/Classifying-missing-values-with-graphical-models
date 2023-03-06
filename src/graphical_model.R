#############################################################
#                  Fit the graphical models
#############################################################

# Get required packages & data
library(tidyverse)
source("src/graphical_model_utils.R")
load("data/multi_normal.RData")
load("data/mixed_data.RData")


fitted_graph <- apply_graph_fitting(
  multi_normal_data_list,
  experiment = "mnar_x2",
  rel_missingness = 0.3,
  alpha = 0.01,
  type = "mvpc"
)
fitted_graph <- apply_graph_fitting(
  mixed_data_list,
  experiment = "mar",
  rel_missingness = 0.6,
  alpha = 0.01,
  type = "pc"
)

# fit and visualize all at the same time
multi_normal_data_list <- fit_and_visualize_graphs(
  multi_normal_data_list,
  type = "pc",
  plot_graph = FALSE
)

# detect missingness for the full list
multi_normal_data_list <- add_classification(multi_normal_data_list)


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
    classification_results$correct_classification[index] <- mean(
      multi_normal_data_list[[experiment]]$detected_missingness_type[[rel_missing_i]] ==
        multi_normal_data_list[[experiment]]$missing_type
    )
    index <- index + 1
  }
}

classification_results %>%
  mutate(
    relative_missingness = as.factor(relative_missingness),
    correct_classification = round(correct_classification, 2),
    experiment = str_to_upper(str_replace(experiment, "_", " "))
  ) %>%
  ggplot(
    aes(x = relative_missingness, y = experiment, fill = correct_classification, label = correct_classification)
  ) +
  geom_tile(col = "white") +
  geom_text(col = "white", size = 5) +
  # scale_fill_manual(
  #   values = c("FALSE" = "orange", "TRUE" = "darkgreen"),
  #   name = "Classification",
  #   labels = c("TRUE" = "Correct", "FALSE" = "Wrong")
  # ) +
  labs(x = "Relative Missingness", y = "Experiment", fill = "Correctly Classified") +
  theme_minimal()


#How the graphs should look like
# Graph without missingness
complete_data <- as.data.frame(mixed_data_list$raw[[1]])[,-1]
suffStat <- list(C = cor(complete_data), n = nrow(complete_data))
pc_algorithm <- pc(suffStat, indepTest = gaussCItest,
                   p = ncol(complete_data), alpha = 0.01)
plot(pc_algorithm, main = "")
#MAR
library(ggm)
graph_mar <- DAG(mixed_data_list$mar$data[[1]]$y ~ 
                   mixed_data_list$mar$data[[1]]$x1
                 + mixed_data_list$mar$data[[1]]$x3
                 + mixed_data_list$mar$data[[1]]$x4 
                 + mixed_data_list$mar$data[[1]]$x5
                 + mixed_data_list$mar$data[[1]]$x6 
                 + mixed_data_list$mar$data[[1]]$x7
                 + mixed_data_list$mar$data[[1]]$x8 
                 + mixed_data_list$mar$data[[1]]$x9,
                   mixed_data_list$mar$data[[1]]$x3 ~
                   mixed_data_list$mar$data[[1]]$x2,
                   mixed_data_list$mar$data[[1]]$x4 ~
                   mixed_data_list$mar$data[[1]]$x3,
                   mixed_data_list$mar$data[[1]]$x8 ~
                   mixed_data_list$mar$data[[1]]$x7,
                   mixed_data_list$mar$data[[1]]$missing_x2 ~ 
                   mixed_data_list$mar$data[[1]]$x1,
                   mixed_data_list$mar$data[[1]]$missing_x3 ~ 
                   mixed_data_list$mar$data[[1]]$x1,
                   mixed_data_list$mar$data[[1]]$missing_y~
                   mixed_data_list$mar$data[[1]]$x1,
                   mixed_data_list$mar$data[[1]]$missing_x5 ~ 
                   mixed_data_list$mar$data[[1]]$x1, 
                   mixed_data_list$mar$data[[1]]$missing_y ~
                   mixed_data_list$mar$data[[1]]$x4,
                   mixed_data_list$mar$data[[1]]$missing_x6 ~ 
                   mixed_data_list$mar$data[[1]]$x7,
                   mixed_data_list$mar$data[[1]]$missing_x8 ~ 
                   mixed_data_list$mar$data[[1]]$x7)

graph_mar_object <- as(graph_mar, "graphNEL")
graph::nodes(graph_mar_object) <- str_split(as.character(graph::nodes(graph_mar_object)),
                                                "\\$", simplify = T)[,4]
edges_mar <- Rgraphviz::buildEdgeList(graph_mar_object)
igraph_obj_mar <- igraph::graph_from_graphnel(graph_mar_object)
igraph::plot.igraph(
  igraph_obj_mar,
  mark.col = "blue",
  vertex.color = "lightblue",
  vertex.size = 20,
  edge.arrow.size = 0.5,
  label.cex = 12,
)
#Rgraphviz::plot(graph_mar_object)

#MNAR
graph_mnar <- DAG(mixed_data_list$mnar$data[[1]]$y ~ 
                    mixed_data_list$mnar$data[[1]]$x1
                  + mixed_data_list$mnar$data[[1]]$x3
                  + mixed_data_list$mnar$data[[1]]$x4 
                  + mixed_data_list$mnar$data[[1]]$x5
                  + mixed_data_list$mnar$data[[1]]$x6 
                  + mixed_data_list$mnar$data[[1]]$x7
                  + mixed_data_list$mnar$data[[1]]$x8 
                  + mixed_data_list$mnar$data[[1]]$x9,
                  mixed_data_list$mnar$data[[1]]$x3 ~
                    mixed_data_list$mnar$data[[1]]$x2,
                  mixed_data_list$mnar$data[[1]]$x4 ~
                    mixed_data_list$mnar$data[[1]]$x3,
                  mixed_data_list$mnar$data[[1]]$x8 ~
                    mixed_data_list$mnar$data[[1]]$x7,
                   mixed_data_list$mnar$data[[1]]$missing_x2 ~ 
                   mixed_data_list$mnar$data[[1]]$x1,
                   mixed_data_list$mnar$data[[1]]$missing_x5 ~ 
                   mixed_data_list$mnar$data[[1]]$x1,
                   mixed_data_list$mnar$data[[1]]$missing_y ~ 
                   mixed_data_list$mnar$data[[1]]$x1,
                   mixed_data_list$mnar$data[[1]]$missing_x3 ~ 
                   mixed_data_list$mnar$data[[1]]$x2,
                   mixed_data_list$mnar$data[[1]]$missing_x5 ~ 
                   mixed_data_list$mnar$data[[1]]$x4,
                   mixed_data_list$mnar$data[[1]]$missing_x6 ~ 
                   mixed_data_list$mnar$data[[1]]$x7,
                   mixed_data_list$mnar$data[[1]]$missing_x8 ~ 
                   mixed_data_list$mnar$data[[1]]$x7,
                   mixed_data_list$mnar$data[[1]]$missing_y ~ 
                   mixed_data_list$mnar$data[[1]]$x7)

graph_mnar_object <- as(graph_mnar, "graphNEL")
graph::nodes(graph_mnar_object) <- str_split(as.character(graph::nodes(graph_mnar_object)),
                                            "\\$", simplify = T)[,4]
edges <- Rgraphviz::buildEdgeList(graph_mnar_object)
igraph_obj <- igraph::graph_from_graphnel(graph_mnar_object)
igraph::plot.igraph(
  igraph_obj,
  mark.col = "blue",
  vertex.color = "lightblue",
  vertex.size = 20,
  edge.arrow.size = 0.5,
  label.cex = 12,
)
#Rgraphviz::plot(graph_mar_object)


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
complete_data <- as.data.frame(mixed_data_list$raw[[1]])[,-1]

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



library(devtools)
install_github("OmegaPetrazzini/NAsImpute")






