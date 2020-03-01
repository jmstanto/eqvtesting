# Equivalence testing for alternative procedures:
# Equivalence Testing in the Organizational Sciences
#
# Jeff Stanton, February 26, 2020


#' Run equivalence tests on a bivariate correlation.
#'
#' This procedure works with the raw data for the two variables. No provision
#' has been made for missing data and the two vectors of data must be the
#' same length.
#'
#' @param x A numeric vector with no missing data.
#' @param y A numeric vector with no missing data; same length as x.
#' @param lleq A value between -1 and 1 indicating the lower limit of the equivalence interval.
#' @param uleq A value between -1 and 1 indicating the upper limit of the equivalence interval.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' equivalentCorrelation(x=iris$Sepal.Length, y=iris$Sepal.Width, lleq=-0.25, uleq=0.25)
equivalentCorrelation <- function(x, y, lleq, uleq) {

  testDF <- data.frame(x, y)

  #----------------------
  # Inference by intervals
  corrCI <- cor.test(testDF$x, testDF$y, conf.level=0.90)$conf.int
  ibi <- ( corrCI[2] < uleq ) && ( corrCI[1] > lleq )
  # Results show TRUE if observed cor within interval
  if (ibi) {
    cat("Inference by intervals shows correlation value within the equivalence interval.\n\n")
  } else {
    cat("Inference by intervals shows correlation value outside the equivalence interval.\n\n")
  }

  # TOST
  tostOut <- TOSTER::TOSTr(n=250, r=cor(testDF$x,testDF$y), low_eqbound_r=lleq, high_eqbound_r=uleq, alpha=0.05)
  # Results plot observed cor and CI within equivalence interval
  # print(tostOut)


  #----------------------
  # Bayes estimation
  # Run the BayesFactor procedure for estimating correlations
  # This command creates an MCMC posterior distribution.
  cat("\nPreparing Bayesian estimation test.\n")
  corrOutMCMC <- BayesFactor::correlationBF(testDF$y, testDF$x, rscale = "medium",
                               nullInterval = c(uleq,lleq), posterior = TRUE,
                               iterations=10000)
  plot(corrOutMCMC[,"rho"]) # Trace and density plots of posterior estimates
  ut <- length(which(corrOutMCMC[,"rho"] > uleq)) # Estimates above equiv zone
  lt <- length(which(corrOutMCMC[,"rho"] < lleq)) # Estimate below equiv zone
  # max(corrOutMCMC[,"rho"]) # Confirms largest
  # min(corrOutMCMC[,"rho"]) # Confirms smallest
  propEQ <- (1 - (ut + lt)/10000) # Proportion of estimates within equivalence interval

  cat("The proportion of posterior estimates within the equivalence interval is:\n")
  cat(propEQ)
  cat("\n")


  #----------------------
  # Bayes factor
  # This command calculates Bayes factors
  cat("\nPreparing Bayes factor test.\n")
  corrOutBF <- BayesFactor::correlationBF(testDF$y, testDF$x, rscale = "medium",
                             nullInterval = c(lleq,uleq), posterior = FALSE)
  # print(corrOutBF) # Overall output
  nullRatio <- corrOutBF[1]/corrOutBF[2] # Bayes factor in favor of null interval
  print(nullRatio)

  return(invisible(NULL))
}



#' Run equivalence tests on a regression weight.
#'
#' The selection of upper and lower limits of the equivalence band must
#' accord with whether or not the regression was standardized. An unstandardized
#' regression (with B weights) must be accompanied by upper and lower limits that
#' are calibrated to fit the B weight being tested. Choose a standardized regression
#' (e.g., with lm.beta() from the lm.beta package) to produce beta weights so that
#' the upper and lower limits can be calibrated on a scale of -1 to 1.
#'
#' @param lmOut Output object from a run of lm().
#' @param coefname Name of coefficient to examine.
#' @param lleq A value between -1 and 1 indicating the lower limit of the equivalence interval.
#' @param uleq A value between -1 and 1 indicating the upper limit of the equivalence interval.
#' @param beta Use TRUE if the regression is standardized and supply values for lleq and uleq between -1 and 1.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' myDF <- data.frame(y=rnorm(250), x=rnorm(250), z=rnorm(250))
#' out <- lm(y ~ x + z, data=myDF)
#' equivalentWeight(lmOut=out, coefname="x", lleq=-0.20, uleq=0.20)
equivalentWeight <- function(lmOut, coefname, lleq, uleq, beta=TRUE) {

  ibi <- ( confint(lmOut, coefname, level=0.90)[1] > lleq ) &&
    ( confint(lmOut, coefname, level=0.90)[2] < uleq )
  # Test will be TRUE if CI for the slope on X is within
  # the equivalence interval.
  if (ibi) {
    cat("Inference by intervals shows regression weight within the equivalence interval.\n\n")
  } else {
    cat("Inference by intervals shows regression weight value outside the equivalence interval.\n\n")
  }


  #----------------------
  # TOST using one sample t-test on B weight
  sampSize <- length(lmOut$residuals)
  TOSTER::TOSTone.raw(m=lmOut$coefficients[coefname],
              mu=0, sd=summary(lmOut)$coefficients[coefname,"Std. Error"]*sqrt(sampSize),
              n=sampSize, low_eqbound=lleq, high_eqbound=uleq,
              alpha=0.05, plot = TRUE, verbose = TRUE)
  # Plot will show CI for slope on X within the equiv interval.


  #----------------------
  # Bayesian parameter estimation
  testDF <- lmOut$model
  lmOutMCMC <- BayesFactor::lmBF(formula(lmOut), data=testDF, posterior=T, iterations=10000)
  plot(lmOutMCMC[,coefname]) # Posterior density and traceplots of the B weight
  ut <- length(which(lmOutMCMC[,coefname] > uleq)) # Estimates larger than uleq
  lt <- length(which(lmOutMCMC[,coefname] < lleq)) # Estimates smaller than lleq
  #max(lmOutMCMC[,coefname]) # Confirms largest
  #min(lmOutMCMC[,coefname]) # Confirms smallest
  propEQ <- 1 - (ut + lt)/10000 # Proportion of posterior estimates in the null interval
  # Recommend that at least 90% of posterior estimates should fall within the null interval
  cat("The proportion of posterior estimates within the equivalence interval is:\n")
  cat(propEQ)
  cat("\n\n")

  #----------------------
  # Bayes factor
  # Calculate a Bayes factor by residualizing on the other predictors
  # Note that the null interval is here must be expressed in standardized terms
  # because the statistic being tested is the correlation coefficient
  if (beta == TRUE) {
    x <- testDF[,coefname]

    # If the regression contains more than one predictor, residualize on the other predictor(s)
    if (length(lmOut$coefficients) > 2) {
      myForm <- formula(lmOut) # Start with the original formula
      delTerm <- paste(".","~",".","-",coefname) # Assemble the update string
      newModel <- update.formula(myForm, delTerm) # Create a new formula without coefname

      lmOut2 <- lm(newModel, data=testDF)
      y <- residuals(lmOut2)
    } else y <- testDF[,1] # Otherwise, if just one predictor, simply test the correlation between the predictor and outcome

    corrOutBF <- BayesFactor::correlationBF(y, x, rscale = "medium",
                                            nullInterval = c(lleq,uleq), posterior = FALSE)
    #corrOutBF
    print(corrOutBF[1]/corrOutBF[2]) # Bayes factor in favor of null interval
    # Recommend that a Bayes factor of 50 or higher is observed to support the null interval
  } else cat("\nNote that a Bayes factor cannot be calculated for unstandardized values of uleq and lleq.\n")

  return(invisible(NULL))

}


#' Run equivalence tests on an ANOVA weight for a two-level factor.
#' The two-level factor you are testing should be effect coded so that it has a mean of 0 and a standard deviation of 1.
#'
#' @param aovOut Output object from a run of lm().
#' @param coefname Name of coefficient to examine.
#' @param lleq A value between -1 and 1 indicating the lower limit of the equivalence interval.
#' @param uleq A value between -1 and 1 indicating the upper limit of the equivalence interval.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' # Use effect coding on the two predictors to produce standardized
#' # coefficients from the aov/lm model.
#' sampSize = 200
#' testDF <- data.frame(y=rnorm(sampSize),
#'                      x=scale(sample.int(n=2,size=sampSize,replace=TRUE)),
#'                      z=scale(sample.int(n=2,size=sampSize,replace=TRUE)))
#'
#' out <- aov(y ~ z * x, data=testDF) # Run a conventional ANOVA
#' equivalentFactor(out, "z", lleq=-0.20, uleq=0.20)
equivalentFactor <- function(aovOut, coefname, lleq, uleq) {
  if (class(aovOut)[1] != "aov") {
    cat("Please supply an output object from a call to aov.\n")
    return(invisible(NULL))
  }

  testDF <- aovOut$model
  mcoef <- mean(testDF[,coefname])
  sdcoef <- sd(testDF[,coefname])
  if (round(mcoef) != 0 | round(sdcoef) != 1) {
    cat("Use effect coding on .\n")
    return(invisible(NULL))
  }

  lmOut <- lm(formula(aovOut), data=testDF)

  equivalentWeight(lmOut=lmOut, coefname=coefname, lleq=lleq, uleq=uleq, beta=TRUE)

  return(invisible(NULL))
}

