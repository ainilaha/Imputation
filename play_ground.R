


### Implement Matrices
#### Categorical Cross Entropy


categorical_cross_entropy <- function(y_true, pred_prob,epsilon=1e-7){
  loss <- rowMeans(y_true * log(pred_prob + epsilon))
  sum(-loss)
}


categorical_cross_accuracy <- function(y_true, y_imput){
  sum(y_true == y_imput)/length(y_true)
}



#### Binary Cross-Entropy Loss


binary_cross_entropy <- function(y_true, pred_prob, epsilon= 1e-7){
  loss  <- y_true * log(y_pred)  + (1 - y_true) * log(1-pred_prob + epsilon)  
  sum(-loss)
}


#### Mean squared error


mean_square_error <- function(y_true, y_pred){
  mean((y_true - y_pred)^2)
}


### Load Mice

require(mice)
require(lattice)
library(dplyr)

library(rMIDAS)
# set_python_env(python ="/opt/anaconda3/bin/python")
set_python_env(x ="C:\\ProgramData\\Anaconda3\\",type = "conda")


### Load Data

head(fdgs)

#### drop colunms

data <- fdgs %>% select(-c(id,wgt.z,hgt.z)) %>% na.omit()
head(data)


#### Create Missing Data

miss_data <- add_missingness(data, prop = 0.1)
# view miss number of miss data by columns
print(sapply(miss_data, function(x) sum(is.na(x))))
# md.pattern(miss_data)


#### Imputing Data with MICE

miss_data <- as.data.frame(add_missingness(data, prop = 0.1))
imp <-  mice(miss_data, print=F)
meth <- imp$meth
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'norm'

imp <- mice(miss_data, m=10, meth = meth, print=F)
imp20  <-  mice.mids(imp, maxit=1, print=F)
plot(imp20)
impt_mice_data <- mice::complete(imp20)


#### Imputing Data with RMIDAS

col_bin <- c('sex')
col_cat <- c('reg') 
# Apply rMIDAS preprocessing steps
data_conv <- rMIDAS::convert(miss_data, 
                             bin_cols = col_bin, 
                             cat_cols = col_cat,
                             minmax_scale = TRUE)

# Train the model for 20 epochs
rmidas_train <- rMIDAS::train(data_conv,
                              training_epochs = 50,
                              layer_structure = c(128,128),
                              input_drop = 0.75,
                              seed = 89)

# Generate 10 imputed datasets
impt_rmidas_data <- rMIDAS::complete(rmidas_train, m = 10,fast = TRUE)




library(missRanger)
impt_ranger_data <- missRanger(miss_data, num.trees = 100, verbose = 0)





#### plot categorical loss

library(ggplot2)
library(gridExtra)

categorical_entropy_matrix <- function(df,miss_df, impt_df,cols,impt_meth="MICE") 
{
  entory_list <- c()
  accuracy_list <- c()
  
  
  for (col in cols)
  {
    # we only need to compare the missing values
    miss_index <- which(is.na(miss_df[,col]))
    true_col <- as.factor(df[,col])
    impt_col <- as.factor(impt_df[,col])
    

    acc <- categorical_cross_accuracy(as.numeric(true_col[miss_index]),
                                      as.numeric(impt_col[miss_index]))

    accuracy_list <- c(accuracy_list,acc)

    
    # NOTE: the imputed data set may have less  levels than the original data
    true_col <-  model.matrix(~0+true_col)
    impt_col <-  model.matrix(~0+impt_col)
    cross_entropy <- categorical_cross_entropy(true_col[miss_index,],impt_col[miss_index,])
    entory_list <- c(entory_list,cross_entropy)
    result <- paste(col,":","categorical cross entropy=",cross_entropy,
                    "categorical cross accuracy=",acc,
                    "method=",impt_meth)
    print(result)
    
  }
  
  matrix_df <- data.frame(categorical=cols, 
                          cross_entropy=entory_list,
                          cross_accuracy=accuracy_list,
                          methods=c(rep(impt_meth , length(cols)))
  )
  
  matrix_df
}


compare_cross_entropy <- function(data,miss_data,impt_list,cols,methods){
  cat_df <- data.frame(categorical=c(), 
                       cross_entropy=c(),
                       methods=c()
  )
  for (i in 1:length(methods)){
    df <- categorical_entropy_matrix(data,miss_data,impt_list[[i]],cols,methods[i])
    cat_df <- rbind(cat_df,df)
  }
  
  plot1 <- ggplot(cat_df, aes(fill=methods, y=cross_entropy, x=categorical)) +
    geom_bar(position="dodge", stat="identity")+
    xlab("categorical colunms")+ylab("cross entropy loss")
  
  
  plot2 <- ggplot(cat_df, aes(fill=methods, y=cross_accuracy, x=categorical)) +
    geom_bar(position="dodge", stat="identity")+
    xlab("categorical colunms")+ylab("cross accuracy")
  grid.arrange(plot1, plot2, ncol=2)
  
}


impt_ranger_data$reg <- as.integer(impt_ranger_data$reg)
impt_ranger_data$sex <- as.integer(impt_ranger_data$sex)
impt_ranger_data <- as.data.frame(impt_ranger_data)
imputed_dataframes <- list(impt_mice_data,impt_rmidas_data[[1]],impt_ranger_data)
miss_data <- as.data.frame(miss_data)
compare_cross_entropy(data,miss_data,imputed_dataframes,cols=c("sex","reg"),methods = c("MICE","MIDAS","Ranger"))







num_matrix <- function(df,miss_df, impt_df,cols,impt_meth="MICE") 
{
  entory_list <- c()
  for (col in cols)
  {
    # we only need to compare the missing values
    miss_index <- which(is.na(miss_df[,col]))
    mse <- mean_square_error(df[miss_index,col],impt_df[miss_index,col])
    rs <- paste(col,":","Mean Square Error=",mse,"method=",impt_meth)
    print(rs)
  }
  matrix_df <- data.frame(num=cols, mse=mse, methods=c(rep(impt_meth , length(cols))))
  matrix_df
}


comp_num <- function(data,miss_data,impt_list,cols,methods){
  num_df <- data.frame(num=c(), 
                       mse=c(),
                       methods=c()
  )
  for (i in 1:length(methods)){
    df <- num_matrix(data,miss_data,impt_list[[i]],cols,methods[i])
    num_df <- rbind(num_df,df)
  }
  
  ggplot(num_df, aes(fill=methods, y=mse, x=num)) +
    geom_bar(position="dodge", stat="identity")+
    xlab("Numerical Colunms")+ylab("Mean Square Error")
  
}
comp_num(data,miss_data,imputed_dataframes,cols=c( "age", "hgt","wgt"),methods = c("MICE","MIDAS","Ranger"))




for (col in c("sex","reg")){
  print(col)
  miss_data[,col]
  }




# 
# y_true <- rbind(c(0, 1, 0), c(0, 0, 1))
# y_pred <- rbind(c(0.05, 0.95, 0), c(0.1, 0.8, 0.1))
# cat_cross_entropy(y_true, y_pred)
# 
# 
# 
# 
# y_true = c(0, 1, 0, 0)
# y_pred = c(0.6, 0.51, 0.94, 0.8)
# 
# bin_cross_entropy(y_true, y_pred)
# mean_square_error(y_true, y_pred)

