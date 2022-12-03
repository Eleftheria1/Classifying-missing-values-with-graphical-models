#############################################################
#                  Fit the graphical models
#############################################################

# Get required packages & data
library(tidyverse)
load("data/multi_normal.RData")
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
    vertex.size = 50,
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

fitted_mvpc_graph <- apply_mvpc(
  multi_normal_data_list,
  experiment = "mar_x2",
  rel_missingness = 0.3,
  alpha = 0.01
)

# fit and visualize all at the same time
for (experiment in names(multi_normal_data_list)[2:length(multi_normal_data_list)]) {
  for (rel_missing_i in seq_along(multi_normal_data_list[[experiment]]$relative_missingness)) {
    fit_mvpc_graph(
      data = multi_normal_data_list[[experiment]]$data[[rel_missing_i]],
      relative_missingness = multi_normal_data_list[[experiment]]$relative_missingness[[rel_missing_i]],
      experiment = experiment,
      plot_graph = TRUE
    )
  }
}



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


missing_data <- as.data.frame(multi_normal_data_list$mar_x2[[1]][[1]])[,-1]
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

mnar_data <- as.data.frame(multi_normal_data_list$mnar_x2[[1]][[3]])[,-1]
suff_missing_mnar <- list(data = mnar_data)
mvpc_mnar <- mvpc(suffStat = suff_missing_mnar, indepTest =  gaussCItest.td,
                 corrMethod =  gaussCItest.permc,
                 p = 6, alpha = 0.05)
plot(mvpc_mnar)


correlation_mnar <- cor(mnar_data)
cor_complete_cases_only_mnar <- cor(na.omit(mnar_data))
correlation_mnar[2,] <- cor_complete_cases_only_mnar[2,]
correlation_mnar[,2] <- cor_complete_cases_only_mnar[,2]
correlation_mnar[5,] <- cor_complete_cases_only_mnar[5,]
correlation_mnar[,5] <- cor_complete_cases_only_mnar[,5]
suffStat_mnar <- list(C = correlation_mnar, n = nrow(mnar_data))
pc_algorithm_mnar <- pc(suffStat_mnar, indepTest = gaussCItest,
                           p = ncol(mnar_data), alpha = 0.05)


plot(pc_algorithm_mnar, main = "")

mcar_data <- as.data.frame(multi_normal_data_list$mcar_x2[[1]][[1]])[,-1]
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











