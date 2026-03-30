hwlet <- function(M,L,type='mid'){
# Hilbert transform pair of orthogonal wavelet bases
# given by Common factor construction of Selesnick (2001),
# with perfect reconstruction condition
#
# Inputs   M      number of zeros at z=-1
#          L      degree of fractional delay
#          type   if 'mid' the factorisation of the Bezout solution
#		  is obtained with all roots of absolute magnitude less than 1
#                 if 'min' the factorization is given by 'min'-phase solutions
#	          (see Selesnisk (2001))
#                 if 'const' the wavelet does not satisfy perfect reconstruction
#           (see Achard and Gannaz 2024)
#
# Outputs  h, g   scaling filters of length 2*(M+L)
#                 obtained by Common Factor
#
# Achard, Clausel, Gannaz, Roueff (2020)

n <- seq(0,L-1,1)
t <- 1/2
d <- cumprod(c(1, (L-n)*(L-n-t)/(n+1)/(n+1+t)))

if(type=='const'){
  q <- (2*L+1)*2^(-M-2*L)*sqrt(2)
}else{
  s1 <- sapply(seq(0,2*M,1),choose,n=2*M)
  s2 <- conv(d,rev(d))
  s <- conv(s1,s2)

  N <- M+L
  C <- convmtx(s,2*N-1)
  C <- C[seq(2,dim(C)[1],2),]

  b <- rep(0,2*N-1)
  b[N] <- 1
  r <- solve(C,b)
  q <- sfact(r,type)$poly
}

f <- conv(q,choose(M,0:M))

h <- conv(f,d)
g <- conv(f,rev(d))

return(list(g=g,h=h,tau=(h+1i*g)/sqrt(2)))
}

