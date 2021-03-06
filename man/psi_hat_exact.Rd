\name{psi_hat_exact}
\alias{psi_hat_exact}
\title{
discrete Fourier transform of the wavelet
}
\description{
Computes the discrete Fourier transform of the wavelet associated to the given filter using \code{scaling_function}. The length of the Fourier transform is equal to the length of the grid where the wavelet is evaluated.
}
\usage{
psi_hat_exact(filter,J=10)
}
\arguments{
  \item{filter}{wavelet filter as obtained with \code{scaling_filter}.}
  \item{J}{ 2^J corresponds to the size of the grid for the discretisation of the wavelet. The default value is set to 10.	}

}
\value{
\item{psih}{Values of the discrete Fourier transform of the wavelet.}
\item{grid}{Frequencies where the Fourier transform is evaluated.}
}
\references{
G. Fay, E. Moulines, F. Roueff, M. S. Taqqu (2009) 
Estimators of long-memory: Fourier versus wavelets. \emph{Journal of Econometrics}, vol. 151, N. 2, pages 159-177.

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

S. Achard, I Gannaz (2019)
Wavelet-Based and Fourier-Based Multivariate Whittle Estimation: multiwave. \emph{Journal of Statistical Software}, Vol 89, N. 6, pages 1-31. 

}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{DWTexact}}, \code{\link{scaling_filter}}
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
psi_hat_exact(filter,J=6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

