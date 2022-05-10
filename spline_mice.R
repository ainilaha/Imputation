require(mice)
require(lattice)
library(dplyr)
library(rMIDAS)
library(ggplot2)
library(gridExtra)
library("GGally")
library(gdata)
library(splines)


data <- fdgs |> dplyr::select(-c(id,wgt.z,hgt.z)) |> na.omit()

# data <- data[1:1000,]

miss_data <- add_missingness(data, prop = 0.15)
miss_data <- as.data.frame(miss_data)

iter_n <- 15

# mice.impute.gamlss()

#### Imputing Data with MICE spline
mice.impute.b.spline = function(y, ry, x, ridge=0.00001, ...)
{
  # mice.impute.b.spline
  # Regression imputations of y given x, with a b.spline regression
  # line, and with random draws of the residuals around the line.
  # Bayesian bootstrap
  #
  
  # print(x)
  
  # x = cbind(1, as.matrix(x))
  
  xobs = x[ry,]
  yobs = y[ry]
  n1 = sum(ry)
  n0 = sum(!ry)
  
  # do here the Bayesian bootstap Rubin p. 124
  u = runif(n1-1)
  u = diff(c(0,u[order(u)],1))
  s = sample(n1, n1, replace=TRUE, prob=u)
  # print(s)
  dotxobs = xobs[s,]
  dotyobs = yobs[s]
  # print(head(dotxobs))
  # print(head(dotyobs))
  df = data.frame(dotxobs)
  
  # print("df...................................")
  # print(head(df))
  # print(nrow(dotxobs))
  # print(length(dotyobs))
  # print(paste("s=",s,"dotxobs=",dotxobs,"dotyobs=",dotyobs))
  
  # print(paste("s=",s,"dotxobs=",dim(dotxobs),"dotyobs=",length(dotyobs)))
  # bs_lm = lm(target~.,data = df)
  formu = "target ~ "
  i = 1
  for (c in colnames(df))
  {
    if (i != 1){
      formu = paste(formu," + ")
    }
    i = i+1
    if(c %in% c("age","wgt","hgt")){
      knots = quantile(unique(df[,c]),1:8/9)
      knots = paste(knots,collapse = ', ')
      formu = paste(formu,"bs(",c,", knots = c(",knots,"))",sep = "")
    }else{
      formu = paste(formu,c)
    }
  }
  i = 1
  
  df$target = dotyobs

  print(formu)
  bs_lm = lm(as.formula(formu), data=df)
  
  pred = predict(bs_lm,data.frame(x[!ry,  ]),se=T)
  # print(summary(pred))
  # print(head(pred$fit))
  
  return(pred$fit)
}


imp2 <-  mice(miss_data, print=F)
meth <- imp2$meth
# meth[c('sex')] <- "rf"
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'b.spline'

imp20 <- mice(miss_data, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp2, maxit=iter_n, print=F)
impt_mice_spline_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_spline_data <- append(impt_mice_spline_data,list(impt_mice))
}


#### Imputing Data with MICE rf

imp <-  mice(miss_data, print=F)
meth <- imp$meth
# meth[c('sex')] <- "rf"
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'rf'

imp <- mice(miss_data, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp, maxit=20, print=F)
impt_mice_rf_data <- list()

for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_rf_data <- append(impt_mice_rf_data,list(impt_mice))
}




#### Imputing Data with MICE cart

imp3 <-  mice(miss_data, print=F)
meth <- imp3$meth
meth[c('sex','reg')] <- "cart"
# meth[c('sex')] <- "cart"
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








create_compare_data <- function(df,miss_df,impt_df_list,nas,col,m=5,
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
    # print("-----------------------------------------")
    # print(head(df2))
    # print(head(df))
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





# cln = colnames(data)
# paste("target", "~ bs(",cln[1],")+",cln[2],"+",cln[3],sep = "")



