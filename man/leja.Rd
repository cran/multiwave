\name{leja}
\alias{leja}
\title{
Evaluation of the roots for spectral factorization
}
\description{
This program orders the values x_in (supposed to be the roots of a polynomial) in this way that computing the polynomial coefficients by using the function \code{poly} yields numerically accurate results.
}
\usage{
leja(x_in)
}
\arguments{
  \item{x_in}{ Roots of a polynomial
	}
}
\value{
Reordering of \code{x_in}
}
\references{
I.W. Selesnick (2001) Hilbert transform pairs of wavelet bases, \emph{IEEE Signal Processing Letters}, Vol 8, N.6, pp 170-173.
}
\author{
Matlab codes provided by Markus Lang : <lang@dsp.rice.edu>
in 1993,  Rice University.
R code by Achard, Clausel, Gannaz, Roueff (2017).
}
\seealso{
\code{\link{seprts}}, \code{\link{sfact}},\code{\link{hwlet}}
}
\examples{
z = exp(1i*(1:100)*2*pi/100)
p1 = signal::poly(z)
p2 = signal::poly(leja(z))
}

