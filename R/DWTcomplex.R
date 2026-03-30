DWTcomplex <- function(x,filter,real=FALSE){

  ## Computes the the discrete wavelet decomposition of a signal x
  ## for a Common-Factor wavelet using the cascade algorithm
  ## on the grid of dyadic integer $2^{-\text{J}}$
  ##
  ## 	INPUT	x		Vector of signal values to
  ##		filter  	Filter of a complex common-factor wavelet transform
  ##
  ## 	OUTPUT  dwt		Vector of wavelet coeffcients
  ##         	indmaxband 	Vector containing the maximal index of wavelet coefficients at each scale
  ##		Jmax		Maximal scale
  ##
  ##					based on the paper of Fay, Moulines, Roueff & Taqqu (2009)
  ##                                      Achard, Clausel, Gannaz & Roueff (2020)
  ##________________________________________________________________________________________________

if(real){
  return(DWTexact(x,filter))
} else {

h <- filter$h
g <- filter$g

n <- length(x)
N <- length(h)
hL <- h
hH <- rev(h*(-1)^(0:(N-1)))
gL <- g
gH <- rev(g*(-1)^(0:(N-1)))

tmp <- multiwave::compute_nj(n,N)
nscale <- tmp$nj
Jmax <- tmp$J
indmaxband <- cumsum(nscale)

# Pyramidal algorithm for combinations of (phih, psih) and (phig, psig)
dwt <- rep(0,indmaxband[(length(indmaxband))])

aoldh <- x
aoldg <- x
indmin <- 1
for(iscale in 1:Jmax){
	acurh <- rep(0,nscale[iscale])
	acurg <- rep(0,nscale[iscale])
 	d <- rep(0,nscale[iscale])
	offset <- 1
	for(p in 1:nscale[iscale]){
		acurh[p] <- sum(hL*aoldh[offset:(offset+N-1)])
        acurg[p] <- sum(gL*aoldg[offset:(offset+N-1)])
		dwt[indmin+p-1] <- (sum(hH*aoldh[offset:(offset+N-1)])+1i*sum(gH*aoldg[offset:(offset+N-1)]))
		offset <- offset+2
  	}
  	aoldh <- acurh
  	aoldg <- acurg
  	indmin <- indmaxband[iscale]+1
}

return(list(dwt=dwt/sqrt(2),indmaxband=indmaxband,Jmax=Jmax))
}

}
