\name{hwlet}
\alias{hwlet}
\title{
Common-Factor wavelet filter coefficients
}
\description{
Provides the Hilbert transform pair of orthogonal wavelet bases given by Common factor construction of Selesnick (2001),
with perfect reconstruction condition
}
\usage{
hwlet(M,L,type)
}
\arguments{
\item{M}{Number of vanishing moments.}
\item{L}{Degree of fractional delay.}
\item{type}{Type of factorization of the common factors. If 'mid' the factorisation of the Bezout solution is obtained with all roots of absolute magnitude less than 1. Three possible values are 'min', 'mid', and 'const'. If 'min' the factorization is given by 'min'-phase solutions (see Selesnisk (2001)). If 'const' the wavelet does not satisfy perfect reconstruction (see Achard and Gannaz 2024).}
}
\value{
\item{h}{Real part of the filter (up to a normalization)}
\item{g}{Imaginary part of the filter (up to a normalization)}
\item{tau}{Common-factor filter, defined by (h+i*g)/sqrt(2).}
}
\references{
I.W. Selesnick (2001) Hilbert transform pairs of wavelet bases, \emph{IEEE Signal Processing Letters}, Vol 8, N.6, pp 170-173.

I.W. Selesnick (2002) The design of approximate Hilbert transform pairs of wavelet bases, \emph{IEEE Transactions on Signal Processing}, Vol 50, N.5, pp 1144-1152.

S. Achard, M. Clausel, I. Gannaz, F. Roueff (2020). New results on approximate Hilbert pairs of wavelet filters with common factor structure. \emph{Applied and Computational Harmonic Analysis}, Vol 49, N.3, pp 1025-1045. 

S. Achard, I. Gannaz (2024). Local Whittle estimation with (quasi-)analytic wavelets. \emph{Journal of Time Series Analysis}, Vol 45, pp 421-443.

}
\author{S. Achard and I. Gannaz
}
\seealso{
\code{\link{DWTcomplex}}
}
\examples{
filter_complex <- hwlet(M=4,L=4,type='const')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

