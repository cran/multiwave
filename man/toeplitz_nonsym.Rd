\name{toeplitz_nonsym}
\alias{toeplitz_nonsym}
\title{
Transform a vector in a non symmetric Toeplitz matrix
}
\description{
Transform a vector in a non symmetric Toeplitz matrix
}
\usage{
toeplitz_nonsym(vec)
}
\arguments{
  \item{vec}{input vector.}
}
\value{
the corresponding matrix.
}
\references{
S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

S. Achard, I Gannaz (2019)
Wavelet-Based and Fourier-Based Multivariate Whittle Estimation: multiwave. \emph{Journal of Statistical Software}, Vol 89, N. 6, pages 1-31. 

}
\author{
S. Achard and I. Gannaz
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{scaling_function}}
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
Htmp <- toeplitz_nonsym(filter)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

