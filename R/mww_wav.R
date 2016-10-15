mww_wav <- function(xwav,index,psih,grid_K,LU=NULL){

## Computes the multivariate Fourier Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	 xwav	Vector of wavelet coefficients
##               index	Vector containing the indexes of xwav 
##			where the coefficients change of scale
##               psih	List containing psih$fct, the Fourier transform of the 
##			wavelet mother at values psih$grid
##               LU	bivariate vector (optional) containing 
##			L, the lowest resolution in wavelet decomposition
##               	U, the maximal resolution in wavelet decomposition
##
##	OUTPUT   d	Long-range parameter estimation
##		 cov	Long-run covariance matrix estimation
##
##                                           		Achard & Gannaz (2014)
##_________________________________________________________________________________


xwav <- as.matrix(xwav)
k <- dim(xwav)[2]

d_univ <- rep(0,k)
for(ll in seq(1,k,1)){
	d_univ[ll] <- optimize(f=function(d){mww_wav_eval(d,xwav=xwav[,ll], index=index,LU=LU)},lower=-10,upper=10)$minimum
}
md <- d_univ
if(k>1){
 md <- nlm(f=function(d){mww_wav_eval(d,xwav=xwav,index=index,LU=LU)},d_univ)$estimate
}

mg <- mww_wav_cov_eval(md,xwav,index,psih,grid_K,LU)

list(d=md,cov=mg)

}
