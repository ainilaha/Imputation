require(mice)
require(lattice)
library(tidyverse)

library(rMIDAS)
# set_python_env(python ="/opt/anaconda3/bin/python")
set_python_env(x ="C:\\ProgramData\\Anaconda3\\",type = "conda")

library(ggplot2)
library(gridExtra)
library("GGally")

library(gdata)



data <- fdgs %>% select(-c(id,wgt.z,hgt.z)) %>% na.omit()

# data$test <- apply( data[ ,c("reg","age")] , 1 , paste , collapse = "-" )

# head(data)


#### Create Missing Data
miss_data <- add_missingness(data, prop = 0.1)
miss_data <- as.data.frame(miss_data)
miss_index <- which(is.na(miss_data[,"reg"]))
# view miss number of miss data by coluemns
print(sapply(miss_data, function(x) sum(is.na(x))))



#### Imputing Data with missRanger
library(missRanger)

impt_ranger_data <- replicate(
  10, 
  as.data.frame(missRanger(miss_data, verbose = 0, num.trees = 100)), 
  simplify = FALSE
)


for (i in 1:10){
  impt_ranger_data[[i]][,"sex"] <- as.integer(impt_ranger_data[[i]][,"sex"])
  impt_ranger_data[[i]][,"reg"] <- as.integer(impt_ranger_data[[i]][,"reg"])
}




#### Imputing Data with MICE

imp <-  mice(miss_data, print=F)
meth <- imp$meth
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'rf'

imp <- mice(miss_data, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp, maxit=15, print=F)
impt_mice_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_data <- append(impt_mice_data,list(impt_mice))
}





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
                              layer_structure = c(128,256,128),
                              input_drop = 0.75,
                              seed = 89)

# Generate 10 imputed datasets
impt_rmidas_data <- rMIDAS::complete(rmidas_train, m = 10,fast = TRUE)





plot_na_pie <- function(col){
  miss_data <- as.data.frame(miss_data)
  miss_index <- which(is.na(miss_data[,col]))
  na_count <- apply(miss_data[miss_index,], 1, function(x) sum(is.na(x)))
  
  temp_df <- miss_data[miss_index,]
  # test_df$na_count<-na_count
  
  na_pattern <- md.pattern(temp_df)
  print(length(miss_index))
  
  na_row_count <- c()
  for (i in levels(as.factor(na_count))){
    cn <- sum(as.numeric(row.names(na_pattern)[which(na_pattern[,ncol(na_pattern)]==i)]))
    na_row_count <- c(na_row_count,cn)
  }
  
  # na_row_count_perct <- round(na_row_count/sum(na_row_count)*100,2)
  pie_labels <- paste0(round(100 * na_row_count/sum(na_row_count), 2), "%")
  
  labels <- paste(levels(as.factor(na_count)),na_row_count,sep = ":")
  
  
  pie(na_row_count, labels = pie_labels,
      main = paste(col,"NAs Counts"),col = rainbow(length(na_row_count)))
  legend("topright", labels, cex = 0.8,
         fill = rainbow(length(na_row_count)))
  labels
}








create_compare_data <- function(df,miss_df,impt_df_list,nas,col,m=10,
                                method="mice", sp_impt="sex"){
  # refer:https://cran.r-project.org/web/packages/gdata/vignettes/mapLevels.pdf
  map <- mapLevels(x=factor(df$sex))
  # we only need to compare the missing values
  
  df$sex <- as.factor(as.numeric(df$sex))
  miss_df <- as.data.frame(miss_df)
  miss_index <- which(is.na(miss_df[,col]))
  
  na_count <- apply(miss_df[miss_index,], 1, function(x) sum(is.na(x)))
  
  df <- df[miss_index,]
  df["source"] <- rep("True",length(miss_index))
  df$na_count <- rep("True(0 na)",length(miss_index))
  for(i in 1:m){
    df2 <- impt_df_list[[i]]
    df2 <- df2[miss_index,]
    df2["source"] <- rep(method,length(miss_index))
    df2$na_count <-  nas[na_count]
    if(sp_impt=="method"){
      df2["source"] <- rep(paste(method,i,sep = "-"),length(miss_index))
    }
    df <- rbind(df2,df)
  }
  # convert integer to boys and girls
  int <- as.integer(df$sex)
  mapLevels(x=int) <- map
  df$sex <- int
  if (sp_impt=="sex"){
    df$source <- apply( df[ ,c("sex","source")] , 1 , paste , collapse = "-" )
  }
  
  # print(head(df))
  df$na_count = as.factor(df$na_count)
  df
}

# 
# df_mice_wgt <- create_compare_data(data,miss_data,impt_mice_data,nas=wgt_nas,
#                                    col = "wgt",method = "mice",
#                                    sp_impt="sex")


