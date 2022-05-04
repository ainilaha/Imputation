require(mice)
require(lattice)
library(tidyverse)
library(rMIDAS)
library(ggplot2)
library(gridExtra)
library("GGally")
library(gdata)
library(ImputeRobust)


data <- fdgs |> dplyr::select(-c(id,wgt.z,hgt.z)) |> na.omit()

# data <- data[1:100,]

miss_data <- add_missingness(data, prop = 0.2)
miss_data <- as.data.frame(miss_data)


library(VIM)
aggr_plot <- aggr(miss_data, col=c('navyblue','red'), 
                  numbers=TRUE, sortVars=TRUE, labels=names(miss_data), 
                  cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

iter_n <- 10

# mice.impute.gamlss()

#### Imputing Data with MICE gamlss

imp2 <-  mice(miss_data, print=F)
meth <- imp2$meth
pred <- imp2$predictorMatrix
pred['sex','age'] <- 0
pred[,'reg'] <- 0
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'gamlssNO'

imp2 <- mice(miss_data, m=10, method = meth, 
             predictorMatrix=pred,
             visitSequence = "monotone",n.cyc = 1)
imp20  <-  mice.mids(imp2, maxit=iter_n)
impt_mice_gamlss_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_gamlss_data <- append(impt_mice_gamlss_data,list(impt_mice))
}




#### Imputing Data with MICE rf

imp <-  mice(miss_data, print=F)
meth <- imp$meth
# meth[c('sex','reg')] <- "rf"
meth[c('sex','reg','age','hgt','wgt')] <- 'rf'

imp <- mice(miss_data, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp, maxit=iter_n, print=F)
impt_mice_rf_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_rf_data <- append(impt_mice_rf_data,list(impt_mice))
}




#### Imputing Data with MICE cart

imp3 <-  mice(miss_data, print=F)
meth <- imp3$meth
meth[c('sex','reg')] <- "cart"
meth[c('age','hgt','wgt')] <- 'cart'

imp3 <- mice(miss_data, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp3, maxit=iter_n, print=F)
impt_mice_cart_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_cart_data <- append(impt_mice_cart_data,list(impt_mice))
}




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
                                method="rf", sp_impt="sex"){
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





imputed.sets <- mice(miss_data, m = 2,
                     method = "gamlss",
                     visitSequence = "monotone",
                     maxit = 1, seed = 973,
                     n.cyc = 1, bf.cyc = 1,
                     cyc = 1)






