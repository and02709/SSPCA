% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiSMD.R
\name{MultiSMD}
\alias{MultiSMD}
\title{This function decomposes a matrix using singular value decomposition ad uses
the SMD function.  This results in vectors corresponding to the loadings
which have undergone shrinkage}
\usage{
MultiSMD(
  x,
  sumabsu,
  sumabsv,
  K = 3,
  niter = 20,
  v,
  trace = TRUE,
  upos,
  uneg,
  vpos,
  vneg
)
}
\arguments{
\item{x}{this is the matrix which is decomposed and the decomposition vectors
corresponding to the loadings are shrunk}

\item{sumabsu}{this parameter determines by how much the decomposition vectors
corresponding to the observations}

\item{sumabsv}{this parameter determines by how much the decomposition vectors
corresponding to the loadings}

\item{K}{determines how many vectors should be calculated}

\item{niter}{number of iterations to allow for convergence}

\item{v}{allows a particular v matrix to be passed to the SMD function}

\item{trace}{option to deisplay potential warnings}

\item{upos}{allows the option to force the vector corresponding to
observations to correspond to the parallel maxima using pmax}

\item{uneg}{allows the option to force the vector corresponding to
observations to correspond to the parallel minima using pmin}

\item{vpos}{allows the option to force the vector corresponding to
loadings to correspond to the parallel maxima using pmax}

\item{vneg}{allows the option to force the vector corresponding to
loadings to correspond to the parallel minima using pmin}
}
\description{
This function decomposes a matrix using singular value decomposition ad uses
the SMD function.  This results in vectors corresponding to the loadings
which have undergone shrinkage
}
\examples{
MultiSMD(x, sumabsu, sumabsv, K=3, niter=20,v, trace=TRUE, upos, uneg, vpos, vneg)
}
\keyword{matrix}
\keyword{multiple}
\keyword{singular}
