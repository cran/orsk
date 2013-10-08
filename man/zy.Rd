\name{zy}
\alias{zy}
\title{
Estimating the Relative Risk Based on the Odds Ratio
}
\description{
Estimating the relative risk based on the (adjusted) Odds Ratio from multiple logistic regression or other multiple regression models. The method was based on Zhang and Yu (JAMA, 1998)
}
\usage{
zy(risk, oddsratio)
}
\arguments{
  \item{risk}{
the risk rate of having a positive outcome in the control or unexposed group
}
  \item{oddsratio}{
odds ratio estimated from multiple logistic regression or other multiple regression models
}
}
\details{
Primarily for the adjusted odds ratio, the estimated relative risk is given by: 

odds ratio / (1 - risk + risk * odds ratio)
}

\value{
the estimated relative risk 
}
\references{
Zhang J, Yu KF (1998).
What's the relative risk? A method of correcting the odds ratio in cohort studies of common outcomes.
\emph{JAMA}, 280(19), 1690-1.
}
\author{
Zhu Wang
}

\examples{
zy(risk=0.18, oddsratio=2.25)
}
\keyword{ odds ratio }
\keyword{ relative risk }
