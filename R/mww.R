mww <- function(x,filter,LU=NULL){

## Computes the multivariate wavelet Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	x	Data (n*k vector)
##		filter	Wavelet filter
##              LU	bivariate vector (optional) containing 
##			L, the lowest resolution in wavelet decomposition
##               	U, the maximal resolution in wavelet decomposition
##
##	OUTPUT  d	Long-range parameter estimation
##		cov	Long-run covariance matrix estimation
##
##						Achard & Gannaz (2014)
##______________________________________________________________________________


x <- as.matrix(x)
k <- dim(x)[2]

d_univ <- rep(0,k)
for(ll in seq(1,k,1)){
	d_univ[ll] <- optimize(f=function(d){mww_eval(d,x=x[,ll],filter=filter,LU=LU)},lower=-10,upper=10)$minimum
}

md <- d_univ
if(k>1){
 md <- nlm(f=function(d){mww_eval(d,x=x,filter=filter,LU=LU)},d_univ)$estimate
}

mg <- mww_cov_eval(md,x,filter,LU)

list(d=md,cov=mg)

}
