\name{scaling_function}
\alias{scaling_function}
\title{
scaling function and the wavelet function
}
\description{
Computes the scaling function and the wavelet function (for compactly supported wavelet) using the cascade algorithm on the grid of dyadic integer \eqn{2^{-J}}
}
\usage{
scaling_function(filter,J)
}
\arguments{
  \item{filter}{wavelet filter as obtained with \code{scaling_filter}.}
  \item{J}{value of the largest scale.	}

}

\value{
\item{phi}{Scaling function.}
\item{psi}{Wavelet function.}
}
\references{
G. Fay, E. Moulines, F. Roueff, M. S. Taqqu (2009) 
Estimators of long-memory: Fourier versus wavelets. \emph{Journal of Econometrics}, vol. 151, N. 2, pages 159-177. 

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

}
\author{
S. Achard and I. Gannaz
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DWTexact}}, \code{\link{scaling_filter}}
}
\note{
This function was rewritten from an original matlab version by Fay et al. (2009)
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
scaling_function(filter,J=6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

