
library(dplyr)

phesant <- function(df) {
  cnt_data <- nrow(df)
  
  # set default data type as unknown
  data_types <- rep('unknown', ncol(df))
  names(data_types) <- colnames(df)
  
  # **********************continuous started*******************************
  num_data <-
    df |> select_if(is.numeric) |> select_if( ~ !is.integer(.x))
  same_ratio <-
    round(1 - sapply(num_data, function(x)
      n_distinct(x)) / cnt_data, 4)
  
  continous_cols <- names(same_ratio[same_ratio < 0.35])
  
  
  
  num_data <- num_data |> select(-all_of(continous_cols))
  distinct_cnt <- sapply(num_data, function(x)
    n_distinct(x))
  ex_continous_cols <- names(distinct_cnt[distinct_cnt > 20])
  continous_cols <- c(continous_cols, ex_continous_cols)
  
  # remove less than 10 participants
  paticipants_cnt <- sapply(df, function(x)
    sum(!is.na(x)))
  data_types[names(paticipants_cnt[paticipants_cnt < 10])] = "remove"
  
  
  # assign to order and binary
  num_data <- num_data |> select(-all_of(ex_continous_cols))
  distinct_cnt <- sapply(num_data, function(x)
    n_distinct(x))
  bin_cols <- names(distinct_cnt[distinct_cnt <= 2])
  
  order_cols <- names(distinct_cnt[distinct_cnt > 2 & distinct_cnt <= 20])
  
  
  
  # **********************continuous end***********************************
  
  
  # **********************integer started*********************************
  int_data <- df |> select_if(is.integer)
  distinct_cnt <- sapply(int_data, function(x)
    n_distinct(x))
  
  
  bin_cols <- c(bin_cols, names(distinct_cnt[distinct_cnt <= 2]))
  order_cols <- c(order_cols, names(distinct_cnt[distinct_cnt > 2 & distinct_cnt <= 20]))
  
  continous_cols <-
    c(continous_cols, names(distinct_cnt[distinct_cnt > 20]))
  
  # **********************integer end*************************************
  
  
  
  
  # **********************categorical (single) started*******************
  cat_data <-
    df |> select_if(function(x)
      is.factor(x) || is.character(x))
  
  # remove categories than less than 10 participants:
  for (col in colnames(cat_data)) {
    cnt <- cat_data |> count(!!as.name(col))
    cnt <- cnt[cnt$n > 10, ]
    cat_data <- cat_data[cat_data[, col] %in% cnt[, col], ]
  }
  # ordered or un-ordered:
  
  distinct_cnt <- sapply(cat_data, function(x)
    n_distinct(x))
  bin_cols <- c(bin_cols, names(distinct_cnt[distinct_cnt <= 2]))
  
  # NEEDS FURTHER VERY TO BE ORDER
  unorder_cols <- names(distinct_cnt[distinct_cnt > 2])
  
  
  # **********************categorical (single) end************************
  
  data_types[continous_cols] <- 'continuous'
  data_types[bin_cols] <- 'binary'
  data_types[order_cols] <- 'ordered'
  data_types[unorder_cols] <- 'unordered'
  
  data_types
  
}
