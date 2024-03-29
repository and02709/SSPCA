% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ydens.R
\name{ydens}
\alias{ydens}
\title{This function provides visualization for the responses.}
\usage{
ydens(ytrain, ytest, nresp, resp.names)
}
\arguments{
\item{ytrain}{training response dataset}

\item{ytest}{testing response dataset}

\item{nresp}{number of response vectors present in datasets}

\item{resp.names}{gives list of names for each response}

\item{sepAnalysis}{tells whether the eigenvectors should be displayed separately}
}
\description{
This function provides visualization for the responses.
}
\examples{
ydens(ytrain, ytest, nresp, resp.names)
}
\keyword{density}
\keyword{plot}
\keyword{response}
