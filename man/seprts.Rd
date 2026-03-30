\name{seprts}
\alias{seprts}
\title{
Evaluation of the roots for spectral factorization
}
\description{
This program is for spectral factorization.
The roots on the unit circle must have even degree.
Roots with high multiplicity will cause problems,
they should be handled by extracting them prior to
using this program.
}
\usage{
seprts(p, type='mid')
}
\arguments{
  \item{p}{ A polynomial which admits a specral factorization.
}
  \item{type}{ If 'mid' the factorisation of the Bezout solution
		is obtained with all roots of absolute magnitude less than 1.
               If 'min' the factorization is given by 'min'-phase solutions
	        (see Selesnisk (2001)).
	}
}
\value{
The roots of the polynomial \code{p}, separated depending of either they are inside the unit circle or on the unit circle. Technical function for the function \code{sfact}.
}
\references{
I.W. Selesnick (2001) Hilbert transform pairs of wavelet bases, \emph{IEEE Signal Processing Letters}, Vol 8, N.6, pp 170-173.

I.W. Selesnick (2002) The design of approximate Hilbert transform pairs of wavelet bases, \emph{IEEE Transactions on Signal Processing}, Vol 50, N.5, pp 1144-1152.
}
\author{
Matlab codes provided by Selesnick (2001). R code by Achard, Clausel, Gannaz, Roueff (2017).
}
\seealso{
\code{\link{sfact}}, \code{\link{leja}},\code{\link{hwlet}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
\keyword{nonparametric}
