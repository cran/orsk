\name{orsk}
\alias{orsk}
\alias{plot.orsk}
\alias{print.orsk}
\alias{summary.orsk}
\title{
Converting Odds Ratio to Relative Risk in Cohort Studies with Partial Data Information
}
\description{
Converting Odds Ratio to Relative Risk in Cohort Studies with Partial Data Information
}
\usage{
orsk(nctr, ntrt, a=NA, al=NA, au=NA, level = 0.95, type="two-sided", 
method = c("grid","optim"), d=1e-4)
\method{plot}{orsk}(x, type=c("RR", "OR"), digits=2, factor=1, amount=NULL, ...)
\method{print}{orsk}(x, ...)
\method{summary}{orsk}(object, nlist=1:5, ...)
}
\arguments{
  \item{nctr}{
sample size of control group from a published study
}
  \item{ntrt}{
sample size of treatment group from a published study
}
  \item{a}{
estimated odds ratio from a published study
}
  \item{al}{
lower bound of confidence interval from a published study
}
  \item{au}{
upper bound of confidence interval from a published study
}
  \item{level}{
level of confidence interval with default 95\% 
}
  \item{method}{
method for converting the odds ratio to the relative risk with default value "grid"
}
  \item{d}{threshold value (delta in the vignette) to filter out solutions if sum of squares > d. Only used with \code{method="grid"}}
 \item{type}{
type of the objective function with default value "two-sided"; or the type of risk to be plotted. For type="RR", distribution of relative risk among scenarios for which the calculated odds ratio and confidence interval coincide with the published values. For type="OR", distribution of risk of the outcome among scenarios for which the calculated odds ratio and confidence interval coincide with the published values.}
 \item{x}{
object of class \code{orsk}
}
 \item{object}{
object of class \code{orsk}
}
 \item{nlist}{ maximum number of solutions displayed }
 \item{digits}{ rounding accuracy for all the numbers given in the published study, with default value 2 }
 \item{factor, amount}{ arguments for scatter plot, see ?jitter function }
  \item{\dots}{ additional arguments for \code{print, summary}. }
}

\details{
Investigators of medical and epidemiological studies are often interested in comparing a risk of a binary outcome 
between a treatment and control group, or between exposed and unexposed. 
Such an outcome can be an onset of a disease or a dichotomized length of labor duration.

From a published study, suppose we are given the information on sample size of control group \code{nctr}, sample size of treatment group \code{ntrt}, estimated odds ratio \code{a}, and confidence interval \code{(al, au)}, how to estimate the relative risk, when the original 2 by 2 contingency table is not directly available? Two methods are proposed to estimate the cells of the contingency table, and to estimate the relative risk. 
}
\value{
An object of class \code{orsk} is returned. The algorithm estimates the number of outcome in control group \code{ctr_yes}, number of outcome free in control group \code{ctr_no}, number of outcome in treatment group \code{trt_yes} and number of outcome free in treatment group \code{trt_no}. Also the results include the corresponding estimated odds ratio with confidence interval, and relative risk and confidence interval, based on the estimated contingency table. 
}
\references{
  Wang, Zhu (2012). Converting Odds Ratio to Relative Risk in Cohort Studies with Partial Data Information. Under review.
 
  Morris, J.A. and Gardner, MJ (1988). Calculating confidence intervals for relative risks (odds ratios) and standardised ratios and rates. \emph{British Medical Journal}, \bold{296}(6632): 1313--1316.

}
\author{
Zhu Wang
}

\examples{
res1 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="grid")
summary(res1)
res2 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au= 3.03, method="optim")
summary(res2)
res3 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, type="lower", method="grid")
summary(res3)
res4 <- orsk(nctr=1636, ntrt=2601, a=2.61, au=3.03, type="upper", method="grid")
summary(res4)
res5 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, type="ci-only", 
method="grid")
summary(res5)
}
\keyword{ odds ratio }
\keyword{ relative risk }
