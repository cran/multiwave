sfact <- function(h,type='mid'){
# spectral factorization of a polynomial h.
#
# Input   h     polynomial
#         type  if 'mid' the factorisation of the Bezout solution
#		is obtained with all roots of absolute magnitude less than 1
#               if 'min' the factorization is given by 'min'-phase solutions
#	        (see Selesnisk (2001))
#
# Outputs b     new polynomial
#         r     roots of new polynomial
#
# # example:
#    g = runif(10)
#    h = conv(g,rev(g))
#    b = sfact(h)$poly
#    h - conv(b,rev(b)) # should be 0

# required subprograms: seprts.m, leja.m
#
# Matlab codes provided by Selesnick
# R code by Achard, Clausel, Gannaz, Roueff (2017)


if(length(h) == 1){
	b <- sqrt(h)
	r <- c()
}else{

# Get the appropriate roots.
r <- seprts(h,type)

# Form the polynomial from the roots
r <- leja(r)
b <- poly(r)
b <- Re(b)


# normalize
b <- b*sqrt(max(h)/sum(abs(b)^2))

if( max(b)+min(b) < 0){
   b <- -b
}

}

return(list(poly = b, roots = r))


}

# ------------------------------------------------------------


seprts <- function(p,type='mid'){
# This program is for spectral factorization.
# The roots on the unit circle must have even degree.
# Roots with high multiplicity will cause problems,
# they should be handled by extracting them prior to
# using this program.
#
# Matlab codes provided by Selesnick
# R code by Achard, Clausel, Gannaz, Roueff (2017)


  SN <- 0.0001    # Small Number (criterion for deciding if a
                  # root is on the unit circle).

  rts <- sort(roots(p),decreasing=TRUE)

  if(type=='min'){
	# The roots INSIDE the unit circle
	r <- rts[abs(rts)<(1-SN)]

	# The roots ON the unit circle
	orts <-  rts[(abs(rts)>=(1-SN)) & (abs(rts)<=(1+SN)) ]
	N <- length(orts)
	if( N>0 ){
	  if( (N %% 2) == 1){
        	cat('Sorry, but there is a problem in seprts.R')
        	r <- c()
	  }else{
	  # Sort roots on the unit circle by angle
	  k <- sort(Arg(orts),index.return=TRUE)$ix
	  orts <- orts[k[seq(1,N,2)]]

	  # Make final list of roots
	  r <- c(r, orts)
	  }
	}
  }else{
  	# The roots INSIDE the unit circle
	irts <- rts[abs(rts)<(1-SN)]

	k <- (abs(Im(irts))<10^-10)
	Iirts <- irts[!k]
	Rirts <- sort(Re(irts[k]))
	Iirts <- sort(Iirts)
	end <- length(Iirts)

	A <- cbind(Iirts[seq(1,end,2)], Iirts[seq(2,end,2)])
	if(end>1){ A1 <- A[seq(1,end/2,2),] }else{ A1 <- c() }
	if(end>3){ A2 <- A[seq(2,end/2,2),] }else{ A2 <- c() }
	Iirts <- c(A2, 1/A1)
	if(length(Rirts)==1){
		Rirts <- 1/Rirts
	}
	if(length(Rirts)>1){
		Rirts <- c(Rirts[seq(2,length(Rirts),2)], 1/Rirts[seq(1,length(Rirts),2)])
	}
	irts <- c(Iirts, Rirts)

	# The roots ON the unit circle
	orts <- rts[(abs(rts)>=(1-SN)) & (abs(rts)<=(1+SN))]
	N <- length(orts)
	if( (N %% 2) == 1 ){
          cat('Sorry, but there is a problem (1) in seprts.R')
          r <- c()
	  return(NULL)
	}

	# The roots at -1
	K<-0
	if(length(orts)>1){
	  f <- find(abs(orts+1)<SN)
	  K <- length(f)
	  if( K > 0 ){
	  	orts[f] <- c()
	  	N <- N - K    # number of roots on unit circle, (besides -1)
	  }
	  if( (K %% 2) == 1 ){
          	cat('Sorry, but there is a problem (2) in seprts.R')
	  	r <- c()
	  	return(NULL)
	  }
	  # The roots at -1 must be separated from the other roots on
	  # the unit circle, otherwise the following sort by angle will
	  # put half the roots at -1 first, and the other half last.

	  # Sort roots on the unit circle by angle
	  angle_sort <- sort(Arg(orts),index.return=TRUE)
	  orts <- orts[angle_sort$ix]

	  # Make final list of roots
	  a <- ( Arg(orts[seq(1,N,2)]) + Arg(orts[seq(2,N,2)]) )/2
	  orts <- exp(1i*a)
	}

	r <- c(irts, orts, rep(-1,K/2))

  }

  return(r)
}


# ------------------------------------------------------------


leja <- function(x_in){
#    Program orders the values x_in (supposed to be the roots of a
#    polynomial) in this way that computing the polynomial coefficients
#    by using the m-file poly yields numerically accurate results.
#    Try, e.g.,
#               z=exp(1i*(1:100)*2*pi/100);
#    and compute
#               p1 = poly(z);
#               p2 = poly(leja(z));
#    which both should lead to the polynomial x^100-1. You will be
#    surprised!
#
#

# Matlab codes provided by Markus Lang : <lang@dsp.rice.edu>
# 1993  Rice University
#
# R code by Achard, Clausel, Gannaz, Roueff (2017)

x <- t(x_in[])
n <- length(x)
x_out <- rep(1,n)

a <- matrix(x[rep(1,n+1),],ncol=n)

#browser()
if(n>0){

  a[1,] <- abs(a[1,])

  ind <- which.max(a[1,1:n])
  dum1 <- a[ind]
  if(ind!=1){
    dum2 <- a[,1]
    a[,1] <- a[,ind]
    a[,ind] <- dum2
  }
  x_out[1] <- a[n,1]

  if( n>1 ){

	a[2,2:n] <- abs(a[2,2:n]-x_out[1])

	if(n>3){
	for(l in 2:n-1){
  		ind <- which.max(prod(a[1:l,l:n]))
  		dum1 <- prod(a[ind])
  		ind <- ind+l-1
  		if(length(ind)!=0){
  		   if(l!=ind){
    			dum2 <- a[,l]
    			a[,l] <- a[,ind]
    			a[,ind] <- dum2
  		}}
		x_out[l] <- a[n,l]
  		a[l+1,(l+1):n] <- abs(a[l+1,(l+1):n]-x_out[l])
	}}
	x_out <- a[n+1,]
	}
	leja <- x_out
}
else{

leja <- c()

}

return(leja)

}



