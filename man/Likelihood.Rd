% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chb_DataAnalysis.R
\name{Likelihood}
\alias{Likelihood}
\title{Return root means square error
Log-likelihood of choosing right}
\usage{
Likelihood(y, pRight, lapseRate = 0)
}
\arguments{
\item{y}{participant choices}

\item{pRight}{probability of choosing right}

\item{lapseRate}{lapse rate}
}
\value{
-2*logLikelihood
I think the adjustment here should multiply all p by
z'(t) = lapseRate + z(t)(1-2*lapseRate)
where z(t) = p_right if subject choice=right
or 1-p_right if subject choice is left
}
\description{
For minimization and by convention, log-likelihood is returned multiplied by -2
}

