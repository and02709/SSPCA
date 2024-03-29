% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QQplot.R
\name{QQplot}
\alias{QQplot}
\title{This function provides a QQ plot for the training and testing data}
\usage{
QQplot(ytrain, ytest, nresp, resp.names)
}
\arguments{
\item{ytrain}{training response dataset}

\item{ytest}{testing response dataset}

\item{nresp}{number of response vectors present in datasets}

\item{resp.names}{gives list of names for each response}
}
\description{
This function provides a QQ plot for the training and testing data
}
\examples{
QQplot(ytrain, ytest, nresp, resp.names)
}
\keyword{QQ}
\keyword{predicted}
\keyword{response}
