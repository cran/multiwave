\name{DWTexact}
\alias{DWTexact}
\title{
Exact discrete wavelet decomposition 
}
\description{
	Computes the discrete wavelet transform of the data using the pyramidal algorithm.
}
\usage{
DWTexact(x, filter)
}
\arguments{
  \item{x}{
 vector of raw data
}
  \item{filter}{Quadrature mirror filter (also called scaling filter, as returned by the \code{scaling_filter} function)
}
}
\value{
\item{dwt}{computable Wavelet coefficients without taking into account the border effect.}
\item{indmaxband}{vector containing the largest index of each
			band, i.e. for \eqn{j > 1}{j > 1} the wavelet coefficients of scale \eqn{j}{j} are \eqn{\code{dwt}(k)}{\code{dwt}(k)} for \eqn{k \in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]}{k in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]} and for  \eqn{j=1}{j=1},  \eqn{\code{dwt}(k)}{\code{dwt}(k)} for \eqn{k \in [1,\code{indmaxband}(1)]}{k in [1,\code{indmaxband}(1)]}.}
\item{Jmax}{largest available scale index (=length of \code{indmaxband}).}
}
\references{
G. Fay, E. Moulines, F. Roueff, M. S. Taqqu (2009) 
Estimators of long-memory: Fourier versus wavelets. \emph{Journal of Econometrics}, vol. 151, N. 2, pages 159-177.

S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.
}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{scaling_filter}}
}
\note{This function was rewritten from an original matlab version by Fay et al. (2009)}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
u <- rnorm(2^10,0,1)
x <- vfracdiff(u,d=0.2)

	resw <- DWTexact(x,filter)
		xwav <- resw$dwt
		index <- resw$indmaxband
		Jmax <- resw$Jmax

## Wavelet scale 1
ws_1 <- xwav[1:index[1]]
## Wavelet scale 2
ws_2 <- xwav[(index[1]+1):index[2]]
## Wavelet scale 3
ws_3 <- xwav[(index[2]+1):index[3]]
### upto Jmax



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
