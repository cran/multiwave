\name{sfact}
\alias{sfact}
\title{
Spectral factorization of a polynomial. 
}
\description{
Spectral factorization of a polynomial h. 
}
\usage{
sfact(h, type='mid')
}
\arguments{
  \item{h}{
polynomial
}
  \item{type}{
	 If 'mid' the factorisation of the Bezout solution
		is obtained with all roots of absolute magnitude less than 1.
               If 'min' the factorization is given by 'min'-phase solutions
	        (see Selesnisk (2001)).
}
}
\value{
\item{poly}{A new polynomial \code{b}, used in the construction of the common-factor filters, such that \code{h - conv(b,rev(b))} is equal to zeros.}
\item{r}{Roots of the polunomial \code{b}}
}
\references{
I.W. Selesnick (2001) Hilbert transform pairs of wavelet bases, \emph{IEEE Signal Processing Letters}, Vol 8, N.6, pp 170-173.

I.W. Selesnick (2002) The design of approximate Hilbert transform pairs of wavelet bases, \emph{IEEE Transactions on Signal Processing}, Vol 50, N.5, pp 1144-1152.

S. Achard, M. Clausel, I. Gannaz, F. Roueff (2020). New results on approximate Hilbert pairs of wavelet filters with common factor structure. \emph{Applied and Computational Harmonic Analysis}, Vol 49, N.3, pp 1025-1045. 

}
\author{
Matlab codes provided by Markus Lang : <lang@dsp.rice.edu>
in 1993,  Rice University.
R code by Achard, Clausel, Gannaz, Roueff (2017).
}
\seealso{
\code{\link{seprts}}, \code{\link{leja}},\code{\link{hwlet}}
}
\examples{
g = runif(10)
h = conv(g,rev(g))
b = sfact(h)$poly
h - conv(b,rev(b)) ## should be zeros
}

