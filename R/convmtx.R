convmtx <- function(v,n){
# Achard, Clausel, Gannaz, Roueff (2017)

mv<-length(v)
nv<-1
if(is.matrix(v)){
	mv <- dim(v)[1]
	nv <- dim(v)[2]
}

v <- as.vector(v)

v0 <- c(v, rep(0,n-1)) 
r <- rep(0,n)
m <- length(v0)
x <- c(r[seq(n,2,-1)], v0)

cidx <- 0:(m-1)
ridx <- seq(n,1,-1)
vt <- array(cidx,dim=c(m,n)) + t(array(ridx,dim=c(n,m))) # Toeplitz subscripts
convol <- array(x[as.vector(vt)],dim=dim(vt)) # actual data

if(mv < nv){
    convol <- t(convol)
}

return(convol)
}
