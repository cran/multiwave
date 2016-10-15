mfw <- function(x,m){

## Computes the multivariate Fourier Whittle estimation for the of 
## the long-memory parameter d and the long-run covariance matrix. 
## 
## 	INPUT	x	Data (n*k matrix)
##		m	Truncation number in Fourier frequencies
##			
##	OUTPUT  d	Long-range parameter estimation
##		cov	Long-run covariance matrix estimation
##
##					based on the paper of Katsumi Shimotsu, 2007
##					Achard & Gannaz (2014)
##__________________________________________________________________________________

x <- as.matrix(x)
k <- dim(x)[2] 

d_univ <- rep(0,k)
for(l in seq(1,k,1)){
	d_univ[l] <- optimize(f=function(d){mfw_eval(d,x=x[,l],m=m)},lower=-10,upper=10)$minimum
}
md <- d_univ
if(k>1){
 md <- nlm(f=function(d){mfw_eval(d,x=x,m=m)},d_univ)$estimate
}

mg <- mfw_cov_eval(md,x,m)

list(d=md,cov=mg)

}
