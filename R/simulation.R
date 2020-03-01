# Equivalence testing for alternative procedures:
# Equivalence Testing in the Organizational Sciences
#
# Jeff Stanton, February 26, 2020


#' Create two random normal populations for use in the simulation.
#'
#' @param popSize A numeric value indicating the size of the population. Defaults to N=10000.
#' @param effectSize A numeric value for the mean difference between the two populations. The procedure creates standard normal populations, so this value is effectively Cohen's d.
#' @return A data frame with two vectors named pop1 and pop2.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' myPops <- setupTwoNormalPops(popSize=50000, effectSize=0.3)
setupTwoNormalPops <- function(popSize=10000, effectSize=0) {

  pop1 <- rnorm(n=popSize,mean=0,   sd=1) # Control population
  pop2 <- rnorm(n=popSize,mean=effectSize, sd=1) # Effect population

  return(data.frame(pop1, pop2))
}


#' Run a simulation for a specified group size using four null/equivalence techniques.
#'
#' @param pop1 A numeric vector of values from which to sample.
#' @param pop2 A numeric vector of values from which to sample.
#' @param cohD A numeric value calibrated as Cohen's d and representing the bandwidth of the null interval.
#' @param groupSize A numeric value indicating how large a group to sample from each population. Defaults to n=200.
#' @param numReps A numeric value indicating how many experiments to run. Defaults to n=2500.
#' @return A list with four elements: A dataframe of all results, a dataframe of TOST details, a data frame of equivalence bounds details, and a data frame of Bayesian estimation details
#' @export
#' @examples
#' myPops <- setupTwoNormalPops(popSize=10000, effectSize=0)
#' outputList <- fourSimulations(myPops$pop1, myPops$pop2, cohD=0.3, groupSize=100, numReps=250)
#' head(outputList[[1]])
fourSimulations <- function(pop1, pop2, cohD=0.3, groupSize=200, numReps=2500) {

  tostOutDF <- NULL
  eqBoundsDF <- NULL
  bfList <- NULL
  hdiBoundsDF <- NULL

  # Do replications of a t.test comparing samples drawn from the two groups
  for (j in 1:numReps) {

    x=sample(pop1,size=groupSize,replace=TRUE)
    y=sample(pop2,size=groupSize,replace=TRUE)

    # Note that TOST provides the conventional t-test results
    tostDF <- data.frame(y=c(x,y),
                         grp=c(rep_len(0,length.out=groupSize), rep_len(1,length.out=groupSize)))
    tostOut <- TOSTER::dataTOSTtwo(tostDF, "y", "grp",
                           var_equal = F, low_eqbound = -cohD, high_eqbound = cohD, desc=F, plots=F)
    tostOutDF <- rbind(tostOutDF, as.data.frame(tostOut$tost))

    # Inference by intervals: Save the data needed
    eqBoundsDF <- rbind(eqBoundsDF, as.data.frame(tostOut$eqb))


    # tostOut$eqb is used to provide upper and lower limits
    # calculated in the metric of the raw data (i.e., taking
    # the cohD effect size and transforming it back into the raw metric
    # using pooled SD from the samples).
    # The alternative which could work ok for standard normal variables
    # would be just to used the cohD as the boundaries.
    tbfOut <- BayesFactor::ttestBF(x, y,
                      mu = 0, nullInterval = c(tostOut$eqb$asDF[,8], tostOut$eqb$asDF[,9]),
                      paired = FALSE, posterior = FALSE)
    bfList <- c(bfList, tbfOut)

    # Now get a posterior sample
    chains <- BayesFactor::ttestBF(x, y, paired = FALSE, posterior = TRUE, iterations=1000, progress=FALSE)
    hdiBoundsDF <- rbind(hdiBoundsDF, as.data.frame(t(as.matrix(quantile(chains[,2],probs=c(0.05,0.95))))))

    if (j == round(j,-1)) cat(".") # Give a dot for every 10 runs
  }

  # Add the sample sizes to the dataframes
  tostOutDF$grpSize <- groupSize
  eqBoundsDF$grpSize <- groupSize
  bfDF <- data.frame(grpSize=tostOutDF$grpSize, ratio=NA) # A new DF for Bayes results
  hdiBoundsDF$grpSize <- groupSize

  # Repair the data frames
  tostOutDF <- data.frame(tostOutDF,check.names = T) # Get rid of names with []
  eqBoundsDF <- data.frame(eqBoundsDF,check.names = T) # Get rid of names with []

  # Calculate additional indicators: Pass fail for TOST and inference by intervals
  tostOutDF$TOSTpass <- (tostOutDF$p.1. < 0.05) & (tostOutDF$p.2. < 0.05)
  eqBoundsDF$EQBpass <- (eqBoundsDF$cil.raw. > eqBoundsDF$low.raw.) & (eqBoundsDF$ciu.raw. < eqBoundsDF$high.raw.)
  # table(tostOutDF$TOSTpass, eqBoundsDF$EQBpass) # SHould confirm that results are identical

  # This extracts the Bayes factor from each run
  for (i in 1:length(tostOutDF$grpSize)) {
    bothbfs <- 1/exp(bfList[[i]]@bayesFactor[["bf"]])
    bfDF$ratio[i] <- bothbfs[2]/bothbfs[1]
    # In certain cases where no data points fall outside of
    # the null interval, the probability for the denominator
    # !(ll < d < ul) may be too small to compute. In this case
    # the resulting ratio is Nan. Replace with an arbitrarily
    # large Bayes factor as a place holder.
    if (is.na(bfDF$ratio[i])) bfDF$ratio[i] = 2.728521e+13
  }
  bfDF$TOSTpass <- tostOutDF$TOSTpass
  bfDF$EQBpass <- eqBoundsDF$EQBpass
  bfDF$logRatio <- log(bfDF$ratio)
  bfDF$or50plus <- (bfDF$ratio >= 50)
  bfDF$hdiLow <- hdiBoundsDF$`5%`
  bfDF$hdiHigh <- hdiBoundsDF$`95%`

  #cat("Mean Bayes Factor for this run: ", mean(bfDF$ratio))
  #cat("Min Bayes Factor for this run: ", min(bfDF$ratio))
  #cat("Median Bayes Factor for this run: ", min(bfDF$ratio))


  # Now process the posterior distributions
  hdiBoundsDF$pass <- (hdiBoundsDF$`5%` > eqBoundsDF$low.raw.) & (hdiBoundsDF$`95%` < eqBoundsDF$high.raw.)
  bfDF$hdiPass <- hdiBoundsDF$pass
  allResults <- bfDF

  return(list(allResults, tostOutDF, eqBoundsDF, hdiBoundsDF))
}


#' Provide a basic diagnostic plot summarizing a run of the simulation.
#'
#' The plot shows the HDI from Bayesian estimation in blue and red, with black
#' lines showing the upper and lower extents of Cohen's d used in the simulation.
#' Any values outside of these lines are non-equivalence according to Bayesian
#' estimation. The log of the Bayes Factor odds ratio is shown on the same graph
#' in green with a green line plotted at log(50). Any green above the line indicates
#' equivalence. Black dots are plotted at the top of the plot each time the
#' TOST and Inference by intervals tests are passed.
#'
#' @param outputList A four item list with a data frame for each item.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' myPops <- setupTwoNormalPops(popSize=10000, effectSize=0)
#' outputList <- fourSimulations(myPops$pop1, myPops$pop2, cohD=0.3, groupSize=100, numReps=250)
#' plotSimulationResults(outputList)
plotSimulationResults <- function(outputList) {
  myDF <- outputList[[1]] # The first data frame in the list is pretty comprehensive

  plot.ts(myDF$hdiLow, col="blue",
          ylim=c(min(myDF$hdiLow),
                 max(myDF$logRatio)),
          xlab="Simulation Run",
          ylab="Cohen's d / Log of Bayes Factor",
          main="Blue/Red:90% HDI; Green:Log Bayes Factor")
  lines(myDF$hdiHigh, col="red")
  lines(myDF$logRatio, col="green")
  abline(h=min(outputList[[3]]$low.cohen.))
  abline(h=min(outputList[[3]]$high.cohen.))
  abline(h=log(50), col="green") # A minimum Bayes Factor of 50 is recommended in the paper

  x = (1:length(myDF$TOSTpass))
  y = myDF$TOSTpass
  dotDF <- data.frame(x,y)
  dotDF <- dotDF[dotDF$y == TRUE,]
  dotDF$y <- max(myDF$logRatio)
  points(dotDF, pch=20)

  return(invisible(NULL))

}



#' Provide basic summary tables summarizing a run of the simulation.
#'
#'
#' @param outputList A four item list with a data frame for each item.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' myPops <- setupTwoNormalPops(popSize=10000, effectSize=0)
#' outputList <- fourSimulations(myPops$pop1, myPops$pop2, cohD=0.3, groupSize=100, numReps=250)
#' tabulateSimulationResults(outputList)
tabulateSimulationResults <- function(outputList) {

  bfDF <- outputList[[1]]
  tostOutDF <- outputList[[2]]
  eqBoundsDF <- outputList[[3]]
  hdiBoundsDF <- outputList[[4]]

  # Compare methods overall
  cat("Summary of TOST results:\n")
  print(table(tostOutDF$TOSTpass))
  cat("\nSummary of Inference-by-Intervals results:\n")
  print(table(eqBoundsDF$EQBpass))
  cat("\nSummary of Bayes Factor results:\n")
  print(table(bfDF$or50plus))
  cat("\nSummary of Bayesian HDI results:\n")
  print(table(hdiBoundsDF$pass))

  cat("\nContingency Tables Comparing Methods:\n")
  cat("\nTOST (row) versus Inference by Intervals (col):\n")
  print(table(tostOutDF$TOSTpass,eqBoundsDF$EQBpass))
  cat("\nTOST (row) versus Bayes Factor (col):\n")
  print(table(tostOutDF$TOSTpass,bfDF$or50plus))
  cat("\nTOST (row) versus Bayesian HDI (col):\n")
  print(table(tostOutDF$TOSTpass,hdiBoundsDF$pass))
  cat("\nBayes Factor (row) versus Bayesian HDI (col):\n")
  print(table(bfDF$or50plus,hdiBoundsDF$pass))

  return(invisible(NULL))
}


