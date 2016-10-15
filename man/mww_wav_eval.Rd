\name{mww_wav_eval}
\alias{mww_wav_eval}

\title{
multivariate wavelet Whittle estimation for data as wavelet coefficients
}
\description{
Evaluates the multivariate wavelet Whittle criterion at a given long-memory parameter vector \code{d} for the already wavelet decomposed data.
}
\usage{
mww_wav_eval(d, xwav, index, LU = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{vector of long-memory parameters (dimension should match dimension of x).
}
  \item{xwav}{
wavelet coefficients matrix (with scales in rows and variables in columns).
}
  \item{index}{
vector containing the largest index of each
			band, i.e. for \eqn{j>1} the wavelet coefficients of scale \eqn{j} are \eqn{\code{dwt}(k)} for \eqn{k \in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]} and for  \eqn{j=1},  \eqn{\code{dwt}(k)} for \eqn{k \in [1,\code{indmaxband}(1)]}.
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
\code{\link{mww}}, \code{\link{mww_cov_eval}},\code{\link{mww_wav}},\code{\link{mww_eval}},\code{\link{mww_wav_cov_eval}}
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

res_mww <- mww_wav_eval(d,xwav,index,LU)
res_d <- optim(rep(0,k),mww_wav_eval,xwav=xwav,index=index,LU=LU,
          method='Nelder-Mead',lower=-Inf,upper=Inf)$par
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
