\name{zy}
\alias{zy}
\title{
Approximate The Relative Risk from the Adjusted Odds Ratio
}
\description{
Estimating relative risk with the adjusted Odds Ratio from multiple logistic regression or other multiple regression models. The method was based on Zhang and Yu (JAMA, 1998)
}
\usage{
zy(risk, oddsratio)
}
\arguments{
  \item{risk}{
the risk rate of having a positive outcome in the control or unexposed group
}
  \item{oddsratio}{
odds ratio obtained from multiple logistic regression or other multiple regression models
}
}
\details{
the approximated relative risk is given by 

odds ratio / (1 - risk + risk * odds ratio)
}

\value{
the approximated relative risk adjusting for the risk rate
}
\references{
Zhang J, Yu KF (1998).
What's the relative risk? A method of correcting the odds ratio in cohort studies of common outcomes.
\emph{JAMA}, \bold{280}(19):1690-1.
}
\author{
Zhu Wang
}

\examples{
zy(risk=0.18, oddsratio=2.25)
}
\keyword{ odds ratio }
\keyword{ relative risk }
