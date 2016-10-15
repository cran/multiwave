fivarma <- function(N, d=0, cov_matrix=diag(length(d)), VAR = NULL, VMA = NULL, skip = 2000){

## Generates N observations of a realisation of a multivariate ARFIMA process X
## according to Model A of Lobato (97).
## Let (e(t)) be a multivariate gaussian process with a covariance matrix cov_matrix.
## The values of the process X are given by the equations:
##		VAR(L).U(t) = VMA(L).e(t)
##	 	diag((1-L)^d).X(t) = U(t)
## where L is the lag-operator.
##
##
## INPUT   N 		Number of time points
##	   d 		Vector of parameters of long-range dependence
##         cov_matrix    Matrix of covariance between the innovations (optional, default is identity)
##	   VAR		Array of VAR coefficient matrices (optional)
##	   VMA		Array of VMA coefficient matrices (optional)
##	   skip		Number of initial observations omitted, after applying the ARMA operator 
##			and the fractional integration (optional, default is 2000)
##         
## OUTPUT  x		Vector containing the N observations of the vector ARFIMA(arlags, d, malags) process
##	   long_run_cov	Matrix of covariance of the spectral density of x around the zero frequency
##			given by Lobato (97)
##
##										Achard & Gannaz (2015)
##__________________________________________________________________________________________________________________


k <- length(d)
if((dim(cov_matrix)[1]!=k)|(dim(cov_matrix)[2]!=k)){ 
	stop('Error, the size covariance matrix does not correspond to the number of long-range parameters.')
}

nn <- N+2*skip


# check if there are AR and MA components
nar <- 0
nma <- 0
if(is.null(VAR)==FALSE){
	nar <- dim(VAR)[3]
	if(is.na(nar)){ 
		nar <-1 
		VAR <- array(VAR,dim=c(k,k,1))
	}
}
if(is.null(VMA)==FALSE){
	nma <- dim(VMA)[3]
	if(is.na(nma)){ 
		nma <-1
		VMA <- array(VMA,dim=c(k,k,1))
	}
}

# tools to compute the long-run covariance matrix
Cvma <- cov_matrix
if(is.null(VMA)==FALSE){
	sumMA <- 0
	for(lag in 1:nma){
		sumMA <- sumMA+VMA[,,lag]
	}
    	Cvma <- (diag(k)+sumMA)%*%Cvma%*%t(diag(k)+sumMA)
}
Cvar <- diag(k)
if(is.null(VAR)==FALSE){
	sumAR <- 0
	for(lag in 1:nar){
		sumAR <- sumAR+VAR[,,lag]
	}
	Cvar <- solve(diag(k)+sumAR)
}


# simulate the data
e <- matrix(rnorm(k*nn),nn)%*%chol(cov_matrix)
long_run_cov <- 0*diag(k)

	u <- varma(N+skip,k,VAR,VMA,cov_matrix,e)
   	long_run_cov <- Cvar%*%Cvma%*%t(Cvar);
	x <- vfracdiff(u,d);
	if(k>1){
	x <- x[(dim(x)[1]-N+1):dim(x)[1],] 
}else{
x <- x[(length(x)[1]-N+1):length(x)[1]] 
}

list(x=x,long_run_cov=long_run_cov)
}

