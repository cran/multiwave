\name{psi_hat_exact_complex}
\alias{psi_hat_exact_complex}
\title{
Discrete Fourier transform of a complex Common-Factor wavelet
}
\description{
Computes the discrete Fourier transform of the complex Common-Factor wavelet associated to the given filters using \code{scaling_function}. The length of the Fourier transform is equal to the length of the grid where the wavelet is evaluated.
}
\usage{
psi_hat_exact_complex(h,g,J=10)
}
\arguments{
\item{h}{Real part of the filter (up to a normalization)}
\item{g}{Imaginary part of the filter (up to a normalization)}
  \item{J}{ 2^J corresponds to the size of the grid for the discretisation of the wavelet. The default value is set to 10.	}

}
\value{
\item{phih}{Values of the discrete Fourier transform of the scaling wavelet.}
\item{psih}{Values of the discrete Fourier transform of the mother wavelet.}
\item{grid}{Frequencies where the Fourier transform is evaluated.}
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
\code{\link{DWTexact}}, \code{\link{scaling_filter}}
}
\examples{
filter_complex <- hwlet(M=4,L=4,type='const')
psi_hat_exact_complex(filter_complex$h,filter_complex$g,J=6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

