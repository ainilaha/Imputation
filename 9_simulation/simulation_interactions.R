

### Data generating functions

makeSurv <- function(n = 2000, loghr = kLogHR) {
  # Creates a survival cohort of n patients. Assumes that censoring is
  # independent of all other variables
  # x1 and x2 are random normal variables
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  # Create the x3 variable
  data$x3 <- 0.5 * (data$x1 + data$x2 - data$x1 * data$x2) + rnorm(n)
  # Underlying log hazard ratio for all variables is the same
  data$y <- with(data, loghr * (x1 + x2 + x3))
  data$survtime <- rexp(n, exp(data$y))
  # Censoring - assume uniform distribution of observation times
  # up to a maximum
  obstime <- runif(nrow(data),
                   min = 0,
                   max = quantile(data$survtime, 0.5))
  data$event <- as.integer(data$survtime <= obstime)
  data$time <- pmin(data$survtime, obstime)
  data$cumhaz <- nelsonaalen(data, time, event)
  # True log hazard and survival time are not seen in the data
  # so remove them
  data$y <- NULL
  data$survtime <- NULL
  return(data)
}


makeMarSurv <- function(data, pmissing = kPmiss){
  # Introduces missing data dependent on event indicator
  # and cumulative hazard and x1 and x2
  logistic <- function(x){
    exp(x) / (1 + exp(x))
  }
  predictions <- function(lp, n){
    # uses the vector of linear predictions (lp) from a logistic model
    # and the expected number of positive responses (n) to generate
    # a set of predictions by modifying the baseline
    trialn <- function(lptrial){
      sum(logistic(lptrial))
    }
    stepsize <- 32
    lptrial <- lp
    while(abs(trialn(lptrial) - n) > 1){
      if (trialn(lptrial) > n){
        # trialn bigger than required
        lptrial <- lptrial - stepsize
      } else {
        lptrial <- lptrial + stepsize
      }
      stepsize <- stepsize / 2
    }
    # Generate predictions from binomial distribution
    as.logical(rbinom(logical(length(lp)), 1, logistic(lptrial)))
  }
  data$x3[predictions(0.1 * data$x1 + 0.1 * data$x2 +
                        0.1 * data$cumhaz + 0.1 * data$event, nrow(data) * pmissing)] <- NA
  return(data)
}
#---------------------------------------

### Functions to analyse data

#---------------------------------------

coxfull <- function(data) {
  # Full data analysis
  coefs <- summary(coxph(myformula, data = data))$coef
  # return a vector of coefficients (est), upper and lower 95% limits
  confint <- cbind(coefs[, 'coef'] - qnorm(0.975) * coefs[, 'se(coef)'], coefs[, 'coef'] + qnorm(0.975) * coefs[, 'se(coef)'])
  out <- cbind(coefs[, 'coef'], confint, kLogHR >= confint[, 1] &
                 kLogHR <= confint[, 2])
  colnames(out) <- c('est', 'lo 95', 'hi 95', 'cover')
  out
}

coximpute <- function(imputed_datasets){
  # Analyses a list of imputed data sets
  docoxmodel <- function(data){
    coxph(myformula, data=data)
  }
  mirafits <- as.mira(lapply(imputed_datasets, docoxmodel))
  out <- summary(pool(mirafits))
  out <- cbind(out, kLogHR >= out[, 'lo 95'] & kLogHR <= out[, 'hi 95'])
  # Whether this confidence interval contains the true hazard ratio
  colnames(out)[length(colnames(out))] <- 'cover'
  out
}

domissf <- function(missdata, reps = 10) {
  # Imputation by missForest
  out <- list()
  for (i in 1:reps) {
    invisible(capture.output(out[[i]] <- missForest(missdata)$ximp))
  }
  out
}

dorfimpute <- function(missdata, reps = 10) {
  # Cumulative hazard is the outcome, filling in missing data
  # using proximities, using a 300-tree Random Forest.
  out <- list()
  for (i in 1:reps) {
    invisible(capture.output(out[[i]] <- rfImpute(cumhaz ~ ., missdata, iter = 10)))
  }
  out
}

mice.impute.rfcont5 <- function(y, ry, x, ...) {
  mice.impute.rfcont(
    y = y,
    ry = ry,
    x = x,
    ntree_cont = 5
  )
}
mice.impute.rfcont10 <- function(y, ry, x, ...) {
  mice.impute.rfcont(
    y = y,
    ry = ry,
    x = x,
    ntree_cont = 10
  )
}
mice.impute.rfcont20 <- function(y, ry, x, ...) {
  mice.impute.rfcont(
    y = y,
    ry = ry,
    x = x,
    ntree_cont = 20
  )
}
mice.impute.rfcont50 <- function(y, ry, x, ...) {
  mice.impute.rfcont(
    y = y,
    ry = ry,
    x = x,
    ntree_cont = 50
  )
}

mice.impute.rfcont100 <- function(y, ry, x, ...) {
  mice.impute.rfcont(
    y = y,
    ry = ry,
    x = x,
    ntree_cont = 100
  )
}
domice <- function(missdata, functions, reps = 10) {
  mids <- mice(
    missdata,
    defaultMethod = functions,
    m = reps,
    visitSequence = 'monotone',
    printFlag = FALSE,
    maxit = 10
  )
  lapply(1:reps, function(x)
    complete(mids, x))
}

doanalysis <- function(x) {
  # Creates a data set, analyses it using different methods, and outputs
  # the result as a matrix of coefficients / SE and coverage
  data <- makeSurv(kSampleSize)
  missdata <- makeMarSurv(data)
  out <- list()
  out$full <- coxfull(data)
  out$rfimpute <- coximpute(dorfimpute(missdata))
  out$missf <- coximpute(domissf(missdata))
  out$rf5 <- coximpute(domice(missdata, c('rfcont5')))
  out$rf10 <- coximpute(domice(missdata, c('rfcont10')))
  out$rf20 <- coximpute(domice(missdata, c('rfcont20')))
  out$rf50 <- coximpute(domice(missdata, c('rfcont50')))
  out$rf100 <- coximpute(domice(missdata, c('rfcont100')))
  out$mice <- coximpute(domice(missdata, c('norm')))
  out
}

#---------------------------------------
### Functions to compare methods
#---------------------------------------

pstar <- function(x) {
  if (x < 0.001) {
    '***'
  } else if (x < 0.01) {
    '**'
  } else if (x < 0.05) {
    '*'
  } else {
    ''
  }
}
compareBias <- function(method1, method2) {
  # Generates a table comparing bias
  # Comparison statistic is the difference in absolute bias
  # (negative means first method is better)
  compareBiasVar <- function(varname) {
    # All coefficients should be kLogHR
    bias1 <- sapply(results, function(x) {
      x[[method1]][varname, 'est']
    }) - kLogHR
    bias2 <- sapply(results, function(x) {
      x[[method2]][varname, 'est']
    }) - kLogHR
    if (sign(mean(bias1)) == -1) {
      bias1 <- -bias1
    }
    if (sign(mean(bias2)) == -1) {
      bias2 <- -bias2
    }
    paste(formatC(
      mean(bias1) - mean(bias2),
      format = 'fg',
      digits = 3
    ),
    pstar(t.test(bias1 - bias2)$p.value))
  }
  sapply(variables, compareBiasVar)
}


compareVariance <- function(method1, method2) {
  # Generates a table comparing precision between two methods
  # Comparison statistic is ratio of variance
  # (smaller means first method is better)
  compareVarianceVar <- function(varname) {
    e1 <- sapply(results, function(x) {
      x[[method1]][varname, 'est']
    })
    e2 <- sapply(results, function(x) {
      x[[method2]][varname, 'est']
    })
    paste(formatC(var(e1) / var(e2), format = 'fg', digits = 3),
          pstar(var.test(e1, e2)$p.value))
  }
  sapply(variables, compareVarianceVar)
}

compareCIlength <- function(method1, method2) {
  # Generates a table comparing coverage percentage between two methods
  # Comparison statistic is the ratio of confidence interval lengths
  # (less than 1 = first better)
  compareCIlengthVar <- function(varname) {
    # Paired t test for bias (difference in estimate)
    len1 <- sapply(results, function(x) {
      x[[method1]][varname, 'hi 95'] -
        x[[method1]][varname, 'lo 95']
    })
    len2 <- sapply(results, function(x) {
      x[[method2]][varname, 'hi 95'] - x[[method2]][varname, 'lo 95']
    })
    paste(formatC(
      mean(len1) / mean(len2),
      format = 'fg',
      digits = 3
    ),
    pstar(t.test(len1 - len2)$p.value))
  }
  sapply(variables, compareCIlengthVar)
}

compareCoverage <- function(method1, method2) {
  # Generates a table comparing coverage percentage between two methods
  # Comparison statistic is the difference in coverage
  # (positive = first better)
  compareCoverageVar <- function(varname) {
    # Paired t test for bias (difference in estimate)
    cov1 <- sapply(results, function(x) {
      x[[method1]][varname, 'cover']
    })
    cov2 <- sapply(results, function(x) {
      x[[method2]][varname, 'cover']
    })
    paste(formatC(100 * (mean(cov1) - mean(cov2)), format = 'f', digits = 1),
          pstar(binom.test(c(
            sum(cov1 == TRUE & cov2 == FALSE),
            sum(cov1 == FALSE & cov2 == TRUE)
          ))$p.value))
  }
  sapply(variables, compareCoverageVar)
}

#---------------------------------------
# Functions to compile and display results
#---------------------------------------

getParams <- function(coef, method) {
  estimates <- sapply(results, function(x) {
    x[[method]][coef, 'est']
  })
  bias <- mean(estimates) - kLogHR
  se_bias <- sd(estimates) / sqrt(length(estimates))
  z <- bias / se_bias
  ci_len <- mean(sapply(results, function(x) {
    x[[method]][coef, 'hi 95'] - x[[method]][coef, 'lo 95']
  }))
  ci_cov <- mean(sapply(results, function(x) {
    x[[method]][coef, 'cover']
  }))
  out <- c(bias, se_bias, z, sd(estimates), ci_len, ci_cov)
  names(out) <- c('bias', 'se_bias', 'z_bias', 'sd', 'ci_len', 'ci_cov')
  out
}

showTable <- function(coef) {
  methods <- c('full',
               'rfimpute',
               'missf',
               'rf5',
               'rf10',
               'rf20',
               'rf50',
               'rf100',
               'mice')
  methodnames <- c(
    'Full data',
    'rfImpute',
    'missForest',
    paste('RF MICE with', c(5, 10, 20, 50, 100), 'trees'),
    'Parametric MICE'
  )
  out <- t(sapply(methods, function(x) {
    getParams(coef, x)
  }))
  out <- formatC(out, digits = 3, format = 'fg')
  out <- rbind(
    c('', 'Standard', 'Z-score ', 'SD of', 'Mean 95%', '95% CI'),
    c(
      'Bias',
      'error of bias',
      'for bias',
      'estimate',
      'CI length',
      'coverage'
    ),
    out
  )
  out <- cbind(c('', '', methodnames), out)
  print(
    xtable(out),
    floating = FALSE,
    include.rownames = FALSE,
    include.colnames = FALSE,
    hline.after = c(0, 2, nrow(out))
  )
}
maketable <- function(comparison) {
  # comparison is a function such as compareCoverage, compareBias
  compare <- cbind(
    comparison('rf10', 'mice'),
    comparison('rf100', 'mice'),
    comparison('rf100', 'rf10')
  )
  compare <- cbind(rownames(compare), compare)
  compare <- rbind(
    c('', 'MICE-RF 10', 'MICE-RF 100', 'MICE-RF 100'),
    c(
      'Coefficient',
      'vs parametric MICE',
      'vs parametric MICE',
      'vs MICE-RF 10'
    ),
    compare
  )
  print(
    xtable(compare),
    include.rownames = FALSE,
    include.colnames = FALSE,
    floating = FALSE,
    hline.after = c(0, 2, nrow(compare))
  )
}