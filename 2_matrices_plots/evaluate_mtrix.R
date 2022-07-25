# library(tensorflow)
# library(keras)
# library(dbplyr)



y_true <- c(2, 1)
y_pred <- rbind(c(0.1, 0.9, 0.8), c(0.05, 0.95, 0))
loss = tf$keras$losses$sparse_categorical_crossentropy(y_true,y_pred)
loss

y_true <- rbind(c(0, 1, 0), c(0, 0, 1))
y_pred <- rbind(c(0.05, 0.95, 0), c(0.1, 0.8, 0.1))
tf$keras$losses$categorical_crossentropy(y_true,y_pred)


# The following method has an overflow/underflow issue on my mackbook pro, it may work on other computers
# library(mltools)
# library(data.table)
# as.matrix(one_hot(as.data.table(data)))


categorical_matrix <- function(df, impt_df,cols) 
{
  # df  : the original data frame
  # impt_df : the imputed data frame
  # cols: categorical columns
  
  # df <- df %>% dplyr::mutate(across(cols,factor))
  # impt_df <- impt_df %>% dplyr::mutate(across(cols,factor))
  # head(impt_df)
  
  for (col in cols)
  {
    true_col <- as.factor(df[,col])
    impt_col <- as.factor(impt_df[,col])
    # NOTE: the imputed data set may have less  levels than the original data
    
    if (0 == length(setdiff(levels(true_col),levels(impt_col))))
    {
      true_col <-  model.matrix(~0+true_col)
      impt_col <-  model.matrix(~0+impt_col)
      # print(class(true_col))
      cross_entropy <- tf$keras$losses$categorical_crossentropy(true_col,impt_col)
      cross_entropy <- sum(as.array(cross_entropy))
      cross_accuracy <- tf$keras$metrics$categorical_accuracy(true_col,impt_col)
      cross_accuracy  <- sum(as.array(cross_accuracy))
      rs <- paste(col,":","categorical crossentropy=",cross_entropy,", categorical accuracy=",cross_accuracy)
      print(rs)
    }
  }
}

categorical_matrix(adult_df,adult_complete[[1]],adult_cat)

numeric_matrix <- function(df, impt_df,cols) 
{
  for (col in  cols)
  {
    mse <- tf$keras$losses$mse(df[,col],impt_df[,col])
    mse <- sum(as.array(mse))
    rs <- paste(col,":","MSE=",mse)
    print(rs)
  }
}

adult_num  <- c("age","fnlwgt", "capital_gain","capital_loss","hours_per_week")
numeric_matrix(adult_df,adult_complete[[1]],adult_num)








df <- adult_df %>% select(adult_cat) %>% mutate(across(adult_cat,factor))
df <- sapply(df, one_hot)

head(df)


wc_true = as.factor(adult_df$native_country[1:5])
wc_miss = as.factor(adult_complete[[1]]$native_country[1:5])
wc_true=model.matrix(~0+wc_true)
wc_miss=model.matrix(~0+wc_miss)
loss = tf$keras$losses$categorical_crossentropy(wc_true,wc_miss)
sum(as.array(loss))

library(cleandata)

adult_df$sex = as.factor(adult_df$sex)
adult_df$class_labels = as.factor(adult_df$class_labels)
adult_cols = encode_binary(adult_df[adult_bin])

# adult_bin_df = adult_df[adult_bin]

adult_bin_df = adult_df %>% select(adult_bin) %>% mutate(across(adult_bin,factor))
adult_bin_misss_df = adult_complete[[1]] %>% select(adult_bin) %>%  mutate(across(adult_bin,factor))

bin_one_hot = model.matrix(~0+adult_bin_misss_df$sex)


tf$keras$losses$binary_crossentropy(adult_bin_df$sex,bin_one_hot)

tf$keras$metrics$binary_accuracy(bin_one_hot,bin_one_hot)


# adult_bin_df = encode_binary(adult_bin_df[adult_bin])


binary_matrix <- function(df, impt_df,cols) 
{

  
  for (col in cols)
  {
    true_col <- as.factor(df[,col])
    impt_col <- as.factor(impt_df[,col])

    # true_col <-  model.matrix(~0+true_col)
    true_col = encode_binary(true_col)
    impt_col <-  model.matrix(~0+impt_col)
    print(class(true_col))
    bin_entropy <- tf$keras$losses$binary_crossentropy(true_col,impt_col)
    bin_entropy <- sum(as.array(bin_entropy))
    binary_accuracy <- tf$keras$metrics$binary_accuracy(true_col,impt_col)
    binary_accuracy  <- sum(as.array(binary_accuracy))
    rs <- paste(col,":","binary_crossentropy=",bin_entropy,", binary_accuracy=",binary_accuracy)
    print(rs)

  }
}

binary_matrix(adult_df,adult_complete[[1]],adult_bin)













miss_cnt = sapply(adult_miss_df, function(x) sum(is.na(x)))



sum(!wc_true==wc_miss)/sum(!adult_miss_df$workclass==wc_true)


