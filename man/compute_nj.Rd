\name{compute_nj}
\alias{compute_nj}
\title{
Wavelets coefficients utilities
}
\description{
Computes the number of wavelet coefficients at each scale.
}
\usage{
compute_nj(n, N)
}
\arguments{
  \item{n}{sample size.}
  \item{N}{filter length.}
}
\value{
\item{nj}{number of coefficients at each scale.}
\item{J}{Number of scales.}
}
\references{
G. Fay, E. Moulines, F. Roueff, M. S. Taqqu (2009) 
Estimators of long-memory: Fourier versus wavelets. \emph{Journal of Econometrics}, vol. 151, N. 2, pages 159-177.

S. Achard, I. Gannaz (2014)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{arXiv}, \code{http://arxiv.org/abs/1412.0391}
}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{DWTexact}}, \code{\link{scaling_filter}}
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
n <- 5^10
N <- length(filter)
compute_nj(n,N)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

