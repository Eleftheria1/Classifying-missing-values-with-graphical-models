library(tidyverse)
load("data/multi_normal.RData")



for (i in 1:length(multi_normal_data_list)) {
    for (k in 1:length(multi_normal_data_list[[i]]$data)) {
        multi_normal_data_list[[i]]$data[[k]]$missing <- if_else(
        is.na(multi_normal_data_list[[i]]$data[[k]]$x2),
        (multi_normal_data_list[[i]]$data[[k]]$missing <- 1),
        (multi_normal_data_list[[i]]$data[[k]]$missing <- 0)
      )
    }
}



