\name{K_eval}
\alias{K_eval}
\title{
Evaluation of function \eqn{K}
}
\description{
Computes the function \eqn{K} as defined in (Achard and Gannaz 2014).
}
\usage{
K_eval(psi_hat,u,d)
}
\arguments{
  \item{psi_hat}{ Fourier transform of the wavelet mother at values \code{u}}

  \item{u}{grid for the approximation of the integral	}

\item{d}{
vector of long-memory parameters. 
}
}
\details{
 \code{K_eval} computes the matrix \eqn{K} with elements
 
\deqn{K(d_l,d_m)=\int u^{(d_l+d_m)} |\code{psi}\_\code{hat}(u)|^2 du }{K(d_l,d_m)=\int u^(d_l+d_m)*|\code{psi}\_\code{hat}(u)|^2 du.}
}
\value{
value of function K as a matrix. 
}
\references{
S. Achard, I. Gannaz (2016)
Multivariate wavelet Whittle estimation in long-range dependence. \emph{Journal of Time Series Analysis}, Vol 37, N. 4, pages 476-512.  \code{http://arxiv.org/abs/1412.0391}.

}
\author{S. Achard and I. Gannaz
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{psi_hat_exact}}
}
\examples{
res_filter <- scaling_filter('Daubechies',8);
filter <- res_filter$h
M <- res_filter$M
alpha <- res_filter$alpha
res_psi <- psi_hat_exact(filter,J=10)
K_eval(res_psi$psih,res_psi$grid,d=c(0.2,0.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

