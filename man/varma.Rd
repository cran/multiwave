\name{varma}
\alias{varma}
\title{
simulation of multivariate ARMA process
}
\description{
generates N observations of a k-vector ARMA process
}
\usage{
varma(N, k = 1, VAR = NULL, VMA = NULL, cov_matrix = diag(k), innov=NULL)
}
\arguments{
  \item{N}{
number of time points.
}
  \item{k}{
dimension of the vectorial ARMA (optional, default is univariate)
}
  \item{VAR}{
array of VAR coefficient matrices (optional).
}
  \item{VMA}{
array of VMA coefficient matrices (optional).
}
 \item{cov_matrix}{
matrix of correlation between the innovations (optional, default is identity).
}
  \item{innov}{
matrix of the innovations (optional, default is a gaussian process).
}
}
\value{
vector containing the N observations of the k-vector ARMA process.
}
\references{
S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.
}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{fivarma}}, \code{\link{vfracdiff}}
}
\examples{

rho1 <- 0.3
rho2 <- 0.8
cov <- matrix(c(1,rho1,rho2,rho1,1,rho1,rho2,rho1,1),3,3)

J <- 9
N <- 2^J
VMA <- diag(c(0.4,0.1,0))
### or another example VAR <- array(c(0.8,0,0,0,0.6,0,0,0,0.2,0,0,0,0,0.4,0,0,0,0.5),dim=c(3,3,2))
VAR <- diag(c(0.8,0.6,0))
x <- varma(N, k=3, cov_matrix=cov, VAR=VAR, VMA=VMA)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{datagen}
