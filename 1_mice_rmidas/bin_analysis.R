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
train_ix = sample(x = 1:nrow(data), size = 2000)
data  = data[train_ix, ]
miss_data <- add_missingness(data, prop = 0.15)
miss_data <- as.data.frame(miss_data)
iter_n <- 15


imp <-  mice(miss_data, print=F)
meth <- imp$meth
meth[c('sex','reg','age','hgt','wgt')] <- 'cart'



imp <- mice(c, m=10, method = meth, print=F)
imp20  <-  mice.mids(imp, maxit=20, print=F)
impt_mice_rf_data <- list()
for (i in 1:10){
  impt_mice <- mice::complete(imp20,action=i)
  impt_mice_rf_data <- append(impt_mice_rf_data,list(impt_mice))
}

one_impt = impt_mice_rf_data[[1]]

miss_index <- which(is.na(miss_data[,"hgt"]))

temp_data = data[miss_index,]
temp_data$err = one_impt[miss_index,"hgt"] - data[miss_index,"hgt"] 

data_g = temp_data |> group_by(age_g = cut(age, breaks=quantile(age, probs=seq(0, 1, by=0.2))))

plot(data_g$age_g,data_g$err)
