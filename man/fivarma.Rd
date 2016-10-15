\name{fivarma}
\alias{fivarma}
\title{
simulation of FIVARMA process
}
\description{
Generates N observations of a realisation of a multivariate FIVARMA process X.
}
\usage{
fivarma(N, d = 0, cov_matrix = diag(length(d)), VAR = NULL,
            VMA = NULL,skip = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
number of time points.
}
  \item{d}{
vector of parameters of long-memory.
}
  \item{cov_matrix}{
matrix of correlation between the innovations (optional, default is identity).
}
  \item{VAR}{
array of VAR coefficient matrices (optional).
}
  \item{VMA}{
array of VMA coefficient matrices (optional).
}
  \item{skip}{
number of initial observations omitted, after applying the ARMA operator and the fractional integration (optional, the default is 2000).
}
}
\details{
Let \eqn{(e(t))_t}{e(t)} be a multivariate gaussian process with a covariance matrix cov_matrix.
The values of the process X are given by the equations:
		\deqn{VAR(L)U(t) = VMA(L)e(t),}{VAR(L)*U(t) = VMA(L)*e(t),} and 


	 	\deqn{diag((1-L)^d)X(t) = U(t)}{diag((1-L)^d)X(t) = U(t)}
		
where \code{L} is the lag-operator.

}
\value{
\item{x}{vector containing the N observations of the vector ARFIMA(arlags, d, malags) process.}
\item{long_run_cov}{matrix of covariance of the spectral density of x around the zero frequency.}
\item{d}{vector of parameters of long-range dependence, modified in case of cointegration.}
}
\references{
R. J. Sela and C. M. Hurvich (2009) Computationaly efficient methods for two multivariate fractionnaly integrated models. \emph{Journal of Time Series Analysis}, Vol 30, N. 6, pages 631-651.

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{varma}}, \code{\link{vfracdiff}}
}
\examples{

rho1 <- 0.3
rho2 <- 0.8
cov <- matrix(c(1,rho1,rho2,rho1,1,rho1,rho2,rho1,1),3,3)
d <- c(0.2,0.3,0.4)

J <- 9
N <- 2^J
VMA <- diag(c(0.4,0.1,0))
### or another example VAR <- array(c(0.8,0,0,0,0.6,0,0,0,0.2,0,0,0,0,0.4,0,0,0,0.5),dim=c(3,3,2))
VAR <- diag(c(0.8,0.6,0))
resp <- fivarma(N, d, cov_matrix=cov, VAR=VAR, VMA=VMA)
x <- resp$x
long_run_cov <- resp$long_run_cov
d <- resp$d



}
\keyword{ts}
\keyword{datagen}
