# eqvtesting
Companion R package for "Evaluating Equivalence and Confirming the Null in the Organizational Sciences" (2021)

Testing and rejecting the null hypothesis is a routine part of quantitative research, but relatively few  researchers prepare for 
*confirming* the null or, similarly, testing a hypothesis of equivalence (e.g., that two group means are practically identical). There are several
statistical testing strategies to accomplish such tests. This R package offers four of those tests and also contains simulator code to generate
the guidance shown in this article:

Stanton, J. M. (2021). Evaluating equivalence and confirming the null in the organizational sciences. Organizational Research Methods, 24(3), 491-512.

The test.R file contains the following test procedures:

equivalentCorrelation() - Tests whether a Pearson's R value is practically equivalent to zero.

equivalentWeight() - Tests whether a B value (slope) from an OLS regression is practically equivalent to zero.

equivalentFactor() - Tests whether a two-level (effect coded) factor in an ANOVA has an effect size that is practically zero.

In addition, the file article.R contains code that generates the tables and figures from the article, including examples of how to use the TOSTER package 
to conduct equivalence tests between the means of two groups.
