library(rMIDAS)
set_python_env(python ="/opt/anaconda3/bin/python")

require(mice)
require(lattice)
library(dplyr)



data <- fdgs %>% select(-c(id,wgt.z,hgt.z)) %>% na.omit()
# head(data)
md.pattern(data)

miss_data <- add_missingness(data, prop = 0.1)
imp <-  mice(miss_data, print=F)
meth <- imp$meth
meth[c('sex','reg')] <- "rf"
meth[c('age','hgt','wgt')] <- 'norm'

imp <- mice(miss_data, m=10, meth = meth, print=F)
imp100  <-  mice.mids(imp, maxit=5, print=F)
plot(imp100)

comp_data <- mice::complete(imp100)





softmax <- function(x){
  exp(x)/sum(exp(x))
  }

eplison_loss <- 1e-7 # to avoid issue raise by log(0):-Inf
cat_cross_entropy <- function(y_true, y_pred){
  loss <- rowMeans(y_true*log(y_pred+eplison_loss))
  sum(-loss)
}


bin_cross_entropy <- function(y_true, y_pred){
  loss  <- y_true* log(y_pred)  + (1-y_true) * log(1-y_pred+eplison_loss)  
  sum(-loss)
}

mean_square_error <- function(y_true, y_pred){
    mean((y_true - y_pred)^2)
}




y_true <- rbind(c(0, 1, 0), c(0, 0, 1))
y_pred <- rbind(c(0.05, 0.95, 0), c(0.1, 0.8, 0.1))
cat_cross_entropy(y_true, y_pred)




y_true = c(0, 1, 0, 0)
y_pred = c(0.6, 0.51, 0.94, 0.8)

bin_cross_entropy(y_true, y_pred)
mean_square_error(y_true, y_pred)

