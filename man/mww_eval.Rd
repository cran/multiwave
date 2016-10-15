\name{mww_eval}
\alias{mww_eval}

\title{
evaluation of multivariate wavelet Whittle estimation
}
\description{
Evaluates the multivariate wavelet Whittle criterion at a given long-memory parameter vector \code{d} using \code{DWTexact} for the wavelet decomposition.
}
\usage{
mww_eval(d, x, filter, LU = NULL)
}
\arguments{
  \item{d}{vector of long-memory parameters (dimension should match dimension of x).
}
  \item{x}{
data (matrix with time in rows and variables in columns).
}
  \item{filter}{
	wavelet filter as obtain with \code{scaling_filter}.
}
  \item{LU}{
	bivariate vector (optional) containing 
			\code{L}, the lowest resolution in wavelet decomposition
               	\code{U}, the maximal resolution in wavelet decomposition.
					(Default values are set to \code{L}=1, and \code{U=Jmax}.)
}
}
\details{
\code{L} is fixing the lower limit of wavelet scales. \code{L} can be increased to avoid finest frequencies that can be corrupted by the presence of high frequency phenomena.

\code{U} is fixing the upper limit of wavelet scales. \code{U} can be decreased when highest frequencies have to be discarded. 
}
\value{
multivariate wavelet Whittle criterion. 
}
\references{
E. Moulines, F. Roueff, M. S. Taqqu (2009) A wavelet whittle estimator of the memory parameter of a nonstationary Gaussian time series. \emph{Annals of statistics}, vol. 36, N. 4, pages 1925-1956

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{mww}}, \code{\link{mww_cov_eval}},\code{\link{mww_wav}},\code{\link{mww_wav_eval}},\code{\link{mww_wav_cov_eval}}
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

## wavelet coefficients definition
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
M <- res_filter$M
alpha <- res_filter$alpha

LU <- c(2,11)

res_mww <- mww_eval(d,x,filter,LU)
k <- length(d)
res_d <- optim(rep(0,k),mww_eval,x=x,filter=filter,
	  	        LU=LU,method='Nelder-Mead',lower=-Inf,upper=Inf)$par

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
