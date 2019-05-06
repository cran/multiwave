\name{mfw}
\alias{mfw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 multivariate Fourier Whittle estimators
}
\description{
Computes the multivariate Fourier Whittle estimators of the long-memory parameters and the long-run covariance matrix also called fractal connectivity. 
}
\usage{
mfw(x, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data (matrix with time in rows and variables in columns).

}
  \item{m}{
truncation number used for the estimation of the periodogram.
}
}
\details{
The choice of m determines the range of frequencies used in the computation of
the periodogram, \eqn{\lambda_j = 2\pi j/N}{lambda_j = 2*pi*j/N}, \eqn{j} = 1,... , \code{m}. The optimal value depends on the spectral properties of the time series such as the presence of short range dependence. In Shimotsu (2007), \code{m} is chosen to be equal to \eqn{N^{0.65}}{N^0.65}.
}
\value{
\item{d}{estimation of the vector of long-memory parameters.}
\item{cov}{estimation of the long-run covariance matrix.}
}
\references{
K. Shimotsu (2007) Gaussian semiparametric estimation of multivariate fractionally integrated processes \emph{Journal of Econometrics} Vol. 137, N. 2, pages 277-310.

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

S. Achard, I Gannaz (2019)
Wavelet-Based and Fourier-Based Multivariate Whittle Estimation: multiwave. \emph{Journal of Statistical Software}, Vol 89, N. 6, pages 1-31. 


}
\author{
S. Achard and I. Gannaz
}

\seealso{
\code{\link{mfw_eval}}, \code{\link{mfw_cov_eval}}
}
\examples{
### Simulation of ARFIMA(0,d,0)
rho <- 0.4
cov <- matrix(c(1,rho,rho,1),2,2)
d <- c(0.4,0.2)
J <- 9
N <- 2^J

resp <- fivarma(N, d, cov_matrix=cov)
x <- resp$x
long_run_cov <- resp$long_run_cov

m <- 57 ## default value of Shimotsu 2007
res_mfw <- mfw(x,m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
