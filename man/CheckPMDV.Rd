% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CheckPMDV.R
\name{CheckPMDV}
\alias{CheckPMDV}
\title{This extracts the desired number of eigenvectors from the X matrix
Performs svd on x matrix in one of two ways.  First is if ncol > nrows,
svd performed on x t(x), otherwise svd performed on t(x) x.  Also, if there
are too many missing values, svd is performed on just x}
\usage{
CheckPMDV(v, x, K)
}
\arguments{
\item{v}{proposed eigenvectors.  Usually set to NULL}

\item{x}{matrix to extract eigenvectors}

\item{K}{number of desired eigenvectors}
}
\description{
This extracts the desired number of eigenvectors from the X matrix
Performs svd on x matrix in one of two ways.  First is if ncol > nrows,
svd performed on x t(x), otherwise svd performed on t(x) x.  Also, if there
are too many missing values, svd is performed on just x
}
\examples{
CheckPMDV(v,x,K)
}
\keyword{decomposition}
\keyword{matrix}
\keyword{penalized}
