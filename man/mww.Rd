\name{mww}
\alias{mww}
\title{
multivariate wavelet Whittle estimation
}
\description{
Computes the multivariate wavelet Whittle estimation for the long-memory parameter  vector \code{d} and the long-run covariance matrix, using \code{DWTexact} for the wavelet decomposition. 
}
\usage{
mww(x, filter, LU = NULL)
}
\arguments{
  \item{x}{
data (matrix with time in rows and variables in columns).
}
  \item{filter}{
	wavelet filter as obtain with \code{scaling_filter}.
}
  \item{LU}{
	bivariate vector (optional) containing 
			\code{L}, the lowest resolution in wavelet decomposition
             	\code{U}, the maximal resolution in wavelet decomposition. (Default values are set to \code{L}=1, and \code{U=Jmax}.)
}
}
\details{
\code{L} is fixing the lower limit of wavelet scales. \code{L} can be increased to avoid finest frequencies that can be corrupted by the presence of high frequency phenomena.

\code{U} is fixing the upper limit of wavelet scales. \code{U} can be decreased when highest frequencies have to be discarded. 
}
\value{
\item{d}{estimation of vector of long-memory parameters.}
\item{cov}{estimation of long-run covariance matrix.}
}
\references{
S. Achard, I. Gannaz (2014)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{arXiv}, \code{http://arxiv.org/abs/1412.0391}
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{mww_eval}}, \code{\link{mww_cov_eval}},\code{\link{mww_wav}},\code{\link{mww_wav_eval}},\code{\link{mww_wav_cov_eval}}
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

## wavelet coefficients definition
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
M <- res_filter$M
alpha <- res_filter$alpha

LU <- c(2,11)

res_mww <- mww(x,filter,LU)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
