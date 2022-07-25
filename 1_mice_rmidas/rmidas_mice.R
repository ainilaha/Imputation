
adult_df = read.csv("https://raw.githubusercontent.com/MIDASverse/MIDASpy/master/Examples/adult_data.csv",
                    # colClasses=c("NULL",NA,NA,NA),
                    row.names = 1
)[1:1500, ]
head(adult_df)

str(adult_df)


library(rMIDAS)
# set_python_env(python ="/opt/anaconda3/bin/python")
set_python_env(x ="C:\\ProgramData\\Anaconda3\\",type = "conda")


set.seed(89)
adult_miss_df <- add_missingness(adult_df, prop = 0.1)
# view miss number of miss data by colunms
sapply(adult_miss_df, function(x) sum(is.na(x)))


adult_cat <- c('workclass','marital_status','relationship','race','education','occupation','native_country')
adult_bin <- c('sex','class_labels')

# Apply rMIDAS preprocessing steps
adult_conv <- rMIDAS::convert(adult_miss_df, 
                      bin_cols = adult_bin, 
                      cat_cols = adult_cat,
                      minmax_scale = TRUE)


# Train the model for 20 epochs
adult_train <- rMIDAS::train(adult_conv,
                     training_epochs = 20,
                     layer_structure = c(128,128),
                     input_drop = 0.75,
                     seed = 89)

# Generate 10 imputed datasets
adult_complete <- rMIDAS::complete(adult_train, m = 10,fast = TRUE)

# Inspect first imputed dataset:
head(adult_complete[[1]])


# Estimate logit model on 10 completed datasets (using Rubin's combination rules)
adult_model <- rMIDAS::combine("class_labels ~ hours_per_week + sex", 
                       adult_complete,
                       family = stats::binomial)

adult_model


require(mice)
require(lattice)

# md.pattern(adult_miss_df)

imp <- mice(adult_miss_df, print=F)
meth <- imp$meth
meth


meth[c('fnlwgt','education_num','capital_gain','capital_loss','hours_per_week')] <- "norm"
meth


### Run imputation
imp <- mice(adult_miss_df, m=5, meth = meth, print=F)
plot(imp)


## Extend Iteration

imp20 <-  mice.mids(imp, maxit=10, print=F)
plot(imp20)

