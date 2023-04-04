#############################################################
#                  Fit the graphical models
#############################################################

# Get required packages & data
library(tidyverse)
source("src/graphical_model_utils.R")
load("data/multi_normal.RData")
load("data/mixed_data.RData")

# VIsualize the individual graphs
fitted_graph <- apply_graph_fitting(
  multi_normal_data_list,
  experiment = "mar_x2",
  rel_missingness = 0.1,
  alpha = 0.01,
  type = "mvpc",
  plot_graph = TRUE
)
fitted_graph <- apply_graph_fitting(
  mixed_data_list,
  experiment = "mcar",
  rel_missingness = 0.1,
  alpha = 0.01,
  type = "mvpc",
  plot_graph = TRUE
)

#-------------------------------------------------------------

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









