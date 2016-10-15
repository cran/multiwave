\name{vfracdiff}
\alias{vfracdiff}
\title{
simulation of vectorial fractional differencing process
}
\description{
Given a vectorial process \code{x} and a vector of long memory parameters \code{d}, this function is producing the corresponding fractional differencing process. 
}
\usage{
vfracdiff(x, d)
}

\arguments{
  \item{x}{
initial process.
}
  \item{d}{
vector of long-memory parameters
}
}
\details{
Given a process \code{x}, this function applied a fractional difference procedure using the formula: 
\deqn{	diag((1-L)^d) x,}{	diag((1-L)^d) x,}

where \code{L} is the lag operator.
}
\value{
vectorial fractional differencing of \code{x}.
}
\references{
S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

K. Shimotsu (2007) Gaussian semiparametric estimation of multivariate fractionally integrated processes \emph{Journal of Econometrics} Vol. 137, N. 2, pages 277-310.
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{varma}}, \code{\link{fivarma}}
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
x <- varma(N, k=3, cov_matrix=cov, VAR=VAR, VMA=VMA)
vx<-vfracdiff(x,d)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{datagen}
