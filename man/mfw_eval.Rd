\name{mfw_eval}
\alias{mfw_eval}

\title{
evaluation of multivariate Fourier Whittle estimator
}
\description{
Evaluates the multivariate Fourier Whittle criterion at a given long-memory parameter value \code{d}.

}
\usage{
mfw_eval(d, x, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
vector of long-memory parameters (dimension should match dimension of x).
}
 \item{x}{data  (matrix with time in rows and variables in columns).

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
multivariate Fourier Whittle estimator computed at point \code{d}. 
}
\references{
K. Shimotsu (2007) Gaussian semiparametric estimation of multivariate fractionally integrated processes \emph{Journal of Econometrics} Vol. 137, N. 2, pages 277-310.

S. Achard, I. Gannaz (2014)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{arXiv}, \code{http://arxiv.org/abs/1412.0391}
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{mfw_cov_eval}}, \code{\link{mfw}}
}
\examples{
### Simulation of ARFIMA(0,d,0)
rho <- 0.4
cov <- matrix(c(1,rho,rho,1),2,2)
d <- c(0.4,0.2)
J <- 9
N <- 2^J

resp <- varfima(N, d, cov_matrix=cov)
x <- resp$x
long_run_cov <- resp$long_run_cov

m <- 57 ## default value of Shimotsu
res_mfw <- mfw(x,m)
d <- res_mfw$d
G <- mfw_eval(d,x,m)
k <- length(d)
res_d <- optim(rep(0,k),mfw_eval,x=x,m=m,method='Nelder-Mead',lower=-Inf,upper=Inf)$par

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
