\name{mww_wav}
\alias{mww_wav}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multivariate wavelet Whittle estimation for data as wavelet coefficients
}
\description{
Computes the multivariate wavelet Whittle estimation of the long-memory parameter  vector \code{d} and the long-run covariance matrix for the already wavelet decomposed data.
}
\usage{
mww_wav(xwav, index, psih, grid_K, LU = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xwav}{
wavelet coefficients matrix (with scales in rows and variables in columns).
}
  \item{index}{
vector containing the largest index of each
			band, i.e. for \eqn{j>1} the wavelet coefficients of scale \eqn{j} are \eqn{\code{dwt}(k)} for \eqn{k \in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]} and for  \eqn{j=1},  \eqn{\code{dwt}(k)} for \eqn{k \in [1,\code{indmaxband}(1)]}.
}
\item{psih}{the Fourier transform of the wavelet mother at values \code{grid_K}.}
\item{grid_K}{the grid for the approximation of the integral in K.}
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
\item{d}{estimation of the vector of long-memory parameters.}
\item{cov}{estimation of the long-run covariance matrix.}
}
\references{
S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

}
\author{
S. Achard and I. Gannaz
}
\seealso{
\code{\link{mww_eval}}, \code{\link{mww_cov_eval}},\code{\link{mww}},\code{\link{mww_wav_eval}},\code{\link{mww_wav_cov_eval}}
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

LU <- c(2,11)

### wavelet decomposition

if(is.matrix(x)){
     N <- dim(x)[1]
     k <- dim(x)[2]
}else{
     N <- length(x)
     k <- 1
}
x <- as.matrix(x,dim=c(N,k))

     ## Wavelet decomposition
     xwav <- matrix(0,N,k)
     for(j in 1:k){
          xx <- x[,j]
             
          resw <- DWTexact(xx,filter)
          xwav_temp <- resw$dwt
          index <- resw$indmaxband
          Jmax <- resw$Jmax
          xwav[1:index[Jmax],j] <- xwav_temp;
     }
     ## we free some memory
     new_xwav <- matrix(0,min(index[Jmax],N),k)
     if(index[Jmax]<N){
          new_xwav[(1:(index[Jmax])),] <- xwav[(1:(index[Jmax])),]
     }
     xwav <- new_xwav
     index <- c(0,index)

##### Compute the wavelet functions 
res_psi <- psi_hat_exact(filter,Jmax)
psih <- res_psi$psih
grid <- res_psi$grid

res_mww <- mww_wav(xwav,index, psih, grid,LU)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
