\name{orsk}
\alias{orsk}
\alias{print.orsk}
\alias{summary.orsk}
\title{
Converting Odds Ratio to Relative Risk with Partial Data Information
}
\description{
Estimating Relative Risk with Odds Ratio and Partial Data Information
}
\usage{
orsk(nctr, ntrt, a=NA, al=NA, au=NA, level = 0.95, type="two-sided", method = c("grid","optim"), d=1e-4)
\method{print}{orsk}(x, ...)
\method{summary}{orsk}(object, nlist=1:5, ...)
}
\arguments{
  \item{nctr}{
total number in control group
}
  \item{ntrt}{
total number in treatment group
}
  \item{a}{
estimated odds ratio
}
  \item{al}{
lower bound of confidence interval
}
  \item{au}{
upper bound of confidence interval
}
  \item{level}{
level of confidence interval with default 95\% 
}
  \item{method}{
estimation method with default "optim"
}
  \item{d}{threshold value to filter out solutions if sum of squares > d. Only used with \code{method="grid"}}
 \item{type}{
type of objective function with the default value "two-sided"}
 \item{x}{
object of class \code{orsk}
}
 \item{object}{
object of class \code{orsk}
}
 \item{nlist}{
an index vector indicating the rank of the smallest sum of squares
}
  \item{\dots}{ additional arguments for \code{print, summary}. }
}

\details{Given nctr, ntrt, odds ratio, and confidence interval, how to obtain relative risk, when n01, n11, n00 and n10 are not directly available? Here we provide two methods to estimate n01 and n11 (thus n00 and n10 can be derived afterwards), and to estimte relative risk from them. The corresponding odds ratio is also compared to the original odds ratio \code{a}.
#             disease (yes)   disease (no)   total
#control         n01              n00        n01+n00=nctr
#treatment       n11              n10        n11+n10=ntrt
}
\value{
An object of class \code{orsk} is returned. The algorithm returns n01, n00, n11, n10, named as ctr_yes, ctr_no, trt_yes, trt_no, respectively. Also the results include the corresponding estimated odds ratio with confidence interval, and relative risk and confidence interval. 
}
\references{
  Lee, S.L., Islam, S., Cassidy, L.D., Abdullah, F., Arca, M.J. (2010). Antibiotics and appendicitis in the pediatric population: An American Pediatric Surgical Association Outcomes and Clinical Trials Committee Systematic Review. \emph{Journal of Pediatric Surgery}, \bold{45}(11):2181--2185. 
 
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
res5 <- orsk(nctr=1636, ntrt=2601, a=2.61, al=2.25, au=3.03, type="ci-only", method="grid")
summary(res5)
}
\keyword{ odds ratio }
\keyword{ relative risk }
