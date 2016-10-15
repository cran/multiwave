\name{scaling_filter}
\alias{scaling_filter}
\title{
wavelet scaling filter coefficients
}
\description{
Computes the filter coefficients of the Haar or Daubechies wavelet family with a specific order
}
\usage{
scaling_filter(family, parameter)
}
\arguments{
  \item{family}{Wavelet family, \code{'Haar'} or \code{'Daubechies'}
}
  \item{parameter}{	Order of the Daubechies wavelet (equal to twice the number of vanishing moments). The value of \code{parameter} can be 2,4,8,10,12,14 and 16.

}
}
\value{
\item{h}{Vector of scaling filter coefficients.}
\item{M}{Number of vanishing moments.}
\item{alpha}{Fourier decay exponent.}
}
\references{
G. Fay, E. Moulines, F. Roueff, M. S. Taqqu (2009) 
Estimators of long-memory: Fourier versus wavelets. \emph{Journal of Econometrics}, vol. 151, N. 2, pages 159-177.

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{DWTexact}}
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
M <- res_filter$M
alpha <- res_filter$alpha
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

