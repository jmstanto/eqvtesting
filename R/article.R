# Article code from:
# Equivalence Testing in the Organizational Sciences
#
# Jeff Stanton, February 26, 2020

#' Show the code from Appendix 1 of the article.
#' The code resides in a Markdown file called appendix1.Rmd that is included
#' with the package. Click on "knit" to create
#' an HTML file with all of the relevant output, or copy/run
#' pieces of code as needed.
#'
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @export
#' @examples
#' showArticleCode()
showArticleCode <- function() {
  codeFile <- system.file("rmd", "appendix1.Rmd", package = "eqvtesting")
  if (rstudioapi::isAvailable()) rstudioapi::navigateToFile(codeFile)
  return(invisible(NULL))
}

#' Display Figure 1 from the article
#'
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @export
#' @examples
#' createFigure1()
createFigure1 <- function() {

  Scenario = NA
  U = NA
  CohensD = NA
  L = NA

  # Data frame with 3 variables: U is upper bound; CohensD is the center, and L is lower bound
  df <- data.frame(Scenario = 1:8,
                   U =       c(0.3, 0.65,  0.65, 0.2,   0.13, -0.1,  0.05, -0.35),
                   CohensD = c(0.0, 0.50,  0.30, 0.15, -0.04, -0.15,-0.30, -0.5),
                   L =      c(-0.3, 0.35, -0.05, 0.1,  -0.21, -0.2, -0.65, -0.65))
  # Plots dots with geom_point and CIs with geom_errorbar
  gOut <- ggplot2::ggplot(df, ggplot2::aes(x = Scenario, y = CohensD)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_errorbar(ggplot2::aes(ymax = U, ymin = L), width=0.25) +
    ggplot2::scale_x_continuous(breaks=1:8) +
    ggplot2::scale_y_continuous(breaks=round(seq(from=-0.6, to=0.6, by=0.1),2))

  print(gOut)
  return(invisible(NULL))
}

#' Show the sequence of demonstration tests from the article
#'
#' @param sigDF The data frame included in the package.
#' @return NULL.
#' @import stats
#' @import graphics
#' @import grDevices
#' @import methods
#' @import utils
#' @export
#' @examples
#' data(sigDF)
#' demonstrateTests(sigDF)
demonstrateTests <- function(sigDF) {

  #-----------------------------------
  # Start by showing some preliminary statistics that are cited
  # or described in the paper
  # Slightly different group sizes
  print(table(sigDF$GENDER))
  # Female   Male
  # 242    298

  # Check the standard deviations by group
  sd(sigDF$SIGPRES[sigDF$GENDER=="Male"])
  # [1] 3.46937
  sd(sigDF$SIGPRES[sigDF$GENDER=="Female"])
  # [1] 3.165619
  # The SDs are slightly different

  sigTout <- t.test(sigDF$SIGPRES[sigDF$GENDER=="Male"], sigDF$SIGPRES[sigDF$GENDER=="Female"])
  print(sigTout) # Not significant
  # Welch Two Sample t-test
  #
  # data:  sigDF$SIGPRES[sigDF$GENDER == "Male"] and sigDF$SIGPRES[sigDF$GENDER == "Female"]
  # t = -0.49734, df = 530.71, p-value = 0.6192
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   -0.7040902  0.4196043
  # sample estimates:
  #   mean of x mean of y
  # 17.81544  17.95768

  # Show the results calibrated as Cohen's d
  cohOut <- effsize::cohen.d(SIGPRES ~ GENDER, data=sigDF)
  print(cohOut)
  print("This confidence interval for d is entirely contained within the +/-0.3 equivalence boundaries.")
  #   Cohen's d
  #
  # d estimate: 0.04262952 (negligible)
  # 95 percent confidence interval:
  #      lower      upper
  # -0.1273733  0.2126323

  # Analyze regular power for this research scenario
  print(pwr::pwr.t.test(n=(min(table(sigDF$GENDER))), d=0.3, type="two.sample"))
  # Two-sample t test power calculation
  #
  # n = 242
  # d = 0.3
  # sig.level = 0.05
  # power = 0.9088086
  # alternative = two.sided
  #
  # NOTE: n is number in *each* group

  # Analyze TOST power for this research scenario
  print(TOSTER::powerTOSTtwo(alpha=0.05, N=min(table(sigDF$GENDER)), low_eqbound_d=-0.3, high_eqbound_d=0.3))
  # The statistical power is 90.21 % for equivalence bounds of -0.3 and 0.3 .
  # [1] 0.902106

  #--------------------------------------------------
  # Demonstrations of equivalence testing start here:

  # Conduct TOST: Gender diffs in SIGPRES
  sigTOST <- TOSTER::dataTOSTtwo(sigDF, "SIGPRES", "GENDER",
                         var_equal = F, low_eqbound = -0.3, high_eqbound = 0.3, desc=T, plots=T)
  print(sigTOST) # Shows that the difference is in fact equivalent to zero

  # Now conduct the ROPE analysis
  # For speed and ease of use, this code uses ttestBF rather than BEST
  print("Preparing Bayesian estimation analysis.")
  bfOut <- BayesFactor::ttestBF(x=sigDF$SIGPRES[sigDF$GENDER=="Male"],
                                y=sigDF$SIGPRES[sigDF$GENDER=="Female"],
                                paired=FALSE,
                                posterior=TRUE,
                                iterations=10000)

  devAskNewPage(ask = FALSE)
  plot(bfOut[,"beta (x - y)"], ask=FALSE)
  print("Lower and upper limits of 90% HDI are:")
  print(quantile(bfOut[,"beta (x - y)"], probs = c(0.05, 0.95)))
  print("The BEST Difference of Means plot shows that at least 95%")
  print("of the posterior estimates fall inside the +/-0.996 equivalence bounds.")

  # The plot shows that plus or minus 0.7 contains 95% of the posterior estimates
  lowCt <- length(which(bfOut[,"beta (x - y)"] < (-0.996))) # How many below the ROPE
  highCt <- length(which(bfOut[,"beta (x - y)"] > (0.996))) # How many above the ROPE
  print("Proportion of posterior estimates within the ROPE:")
  print(1 - (lowCt + highCt)/length(bfOut[,"beta (x - y)"])) # Proportion within the ROPE

  # And finally, the Bayes factor
  tbfOut2 <- BayesFactor::ttestBF(x=sigDF$SIGPRES[sigDF$GENDER=="Male"],
                     y=sigDF$SIGPRES[sigDF$GENDER=="Female"],
                     mu = 0, nullInterval = c(-0.3,0.3),
                     paired = FALSE, posterior = FALSE)
  print("Output from the ttestBF() procedure:")
  print(tbfOut2) # Produces two BF values, one representing the space between d=-0.3 and d=0.3
  # And the other representing all other area. The ratio of these is 2295
  # As an odds ratio it is very strong evidence in favor of the null interval

  # This extracts and divides the two Bayes factors from the output object.
  print("Ratio of two Bayes factors:")
  print(
  exp(tbfOut2@numerator[["Alt., r=0.707 -0.3<d<0.3"]]@analysis[["bf"]])/
    exp(tbfOut2@numerator[["Alt., r=0.707 !(-0.3<d<0.3)"]]@analysis[["bf"]])
  )

  return(invisible(NULL))
}
