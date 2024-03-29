% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSPC.R
\name{SSPC}
\alias{SSPC}
\title{This function supervises the decomposition of the matrix X into vectors
corresponding to the observations and sparse vectors corresponding to
the loadings.  It handles all operations with getting the desired
eigenvectors, and encodes the data.  This information is passed
to subsequent functions}
\usage{
SSPC(
  xtrain,
  ytrain,
  xtest,
  ycont,
  ybinary,
  nresp,
  sumabsv = 4,
  niter = 20,
  K = 1,
  orth = TRUE,
  trace = TRUE,
  v = NULL,
  center = FALSE,
  cnames = NULL,
  vpos = FALSE,
  vneg = FALSE,
  compute.pve = TRUE,
  strictEV = TRUE
)
}
\arguments{
\item{xtrain}{training predictor dataset}

\item{ytrain}{training response dataset}

\item{xtest}{testing predictor dataset}

\item{ycont}{flag for whether the response data is continuous}

\item{ybinary}{flag for whether the response data is binary}

\item{nresp}{number of response vectors in the training response dataset}

\item{sumabsv}{this parameter determines by how much the decomposition vectors
corresponding to the loadings}

\item{niter}{number of iterations to allow for convergence}

\item{K}{determines how many vectors should be calculated}

\item{trace}{option to deisplay potential warnings}

\item{v}{allows a particular v matrix to be passed to PMDL1}

\item{center}{this tells the algorithm whether centering is required
this has been handled in other programs so has been disabled}

\item{cnames}{gives potential column names}

\item{vpos}{allows the option to force the vector corresponding to
loadings to correspond to the parallel maxima using pmax}

\item{vneg}{allows the option to force the vector corresponding to
loadings to correspond to the parallel minima using pmin}

\item{compute.pve}{no longer useful}

\item{strictEV}{determines whether the specified K must be strictly followed
or whether a more appropriate number of vectors can be used}

\item{upos}{allows the option to force the vector corresponding to
observations to correspond to the parallel maxima using pmax}

\item{uneg}{allows the option to force the vector corresponding to
observations to correspond to the parallel minima using pmin}
}
\description{
This function supervises the decomposition of the matrix X into vectors
corresponding to the observations and sparse vectors corresponding to
the loadings.  It handles all operations with getting the desired
eigenvectors, and encodes the data.  This information is passed
to subsequent functions
}
\examples{
SSPC(xtrain, ytrain, xtest, ycont, ybinary, nresp, sumabsv=4, 
  niter=20, K=1, orth=TRUE, trace=TRUE, v=NULL, center=FALSE, cnames=NULL, 
  vpos=FALSE, vneg=FALSE, compute.pve=TRUE, strictEV=TRUE)
}
\keyword{analysis}
\keyword{component}
\keyword{principal}
\keyword{sparse}
\keyword{supervised}
