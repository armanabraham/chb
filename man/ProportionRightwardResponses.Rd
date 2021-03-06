% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chb_DataAnalysis.R
\name{ProportionRightwardResponses}
\alias{ProportionRightwardResponses}
\title{Compute proportion rightward responses}
\usage{
ProportionRightwardResponses(glmData)
}
\arguments{
\item{glmData}{}
}
\value{
Proportion rightward responses and lapse rates for each subject, session (run) and condition
}
\description{
Given single trial data in the form of glmData structure, computes
proportion of choosing the stimulus on the right side (rightward responses).
Also computes lapse rates and merges both rightward responses and lapse rates
together into one dataframe.
}
\examples{
ProportionRightwardResponses(oneSubjectData)

}

