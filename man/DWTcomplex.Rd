\name{DWTcomplex}
\alias{DWTcomplex}
\title{
Exact discrete wavelet decomposition with common-factor wavelets
}
\description{
	Computes the discrete wavelet transform of the data using the pyramidal algorithm.
}
\usage{
DWTcomplex(x, filter, real)
}
\arguments{
  \item{x}{
 vector of raw data
}
  \item{filter}{Common-factor wavelet filters, as returned by the \code{hwlet} function, or real wavelet filters, as returned by the \code{sclaing_filter} function.}
   \item{real}{
Precise if the filter is a real filter (obtained with \code{scaling_filter}) or a complex filter (obtained with \code{hwlet}). The default value is FALSE
}
}
\value{
\item{dwt}{computable Wavelet coefficients without taking into account the border effect.}
\item{indmaxband}{vector containing the largest index of each
			band, i.e. for \eqn{j > 1}{j > 1} the wavelet coefficients of scale \eqn{j}{j} are \eqn{\code{dwt}(k)}{\code{dwt}(k)} for \eqn{k \in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]}{k in [\code{indmaxband}(j-1)+1,\code{indmaxband}(j)]} and for  \eqn{j=1}{j=1},  \eqn{\code{dwt}(k)}{\code{dwt}(k)} for \eqn{k \in [1,\code{indmaxband}(1)]}{k in [1,\code{indmaxband}(1)]}.}
\item{Jmax}{largest available scale index (=length of \code{indmaxband}).}
}
\references{
S. Achard, M. Clausel, I. Gannaz, F. Roueff (2020). New results on approximate Hilbert pairs of wavelet filters with common factor structure. \emph{Applied and Computational Harmonic Analysis}, Vol 49, N.3, pp 1025-1045.

S. Achard, I. Gannaz (2024). Local Whittle estimation with (quasi-)analytic wavelets. \emph{Journal of Time Series Analysis}, Vol 45, pp 421-443.
}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{hwlet}}
}
\examples{
filter_complex <- hwlet(M=4,L=4,type='const')
u <- rnorm(2^10,0,1)
x <- vfracdiff(u,d=0.2)

	resw <- DWTcomplex(x,filter_complex)
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
