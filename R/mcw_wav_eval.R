mcw_wav_eval <- function (d, xwav, index, LU = NULL)
{
    if (is.matrix(xwav)) {
        N <- dim(xwav)[1]
        k <- dim(xwav)[2]
    }
    else {
        N <- length(xwav)
        k <- 1
    }
    Jmax <- length(index) - 1
    if (Jmax<2){
      stop("Unsufficient length of time series")
    }
    if (is.null(LU)) {
        LU <- c(min(3,Jmax-1), Jmax)
    }
    L <- max(LU[1], 1)
    U <- min(LU[2], Jmax)
    nscale <- U - L + 1
    n <- index[U + 1] - index[L]
    sum_xwav <- array(0, dim = c(nscale, k, k))
    vect <- matrix(0, nscale, 1)
    for (f in 1:nscale) {
        j <- L + f - 1
        fj <- 1 - j
        if (k > 1) {
            temp <- xwav[(index[j] + 1):index[j + 1], ] %*% diag(2^(fj * d))
        }
        else {
            temp <- xwav[(index[j] + 1):index[j + 1]] * 2^(fj * d)
        }
        sum_xwav[f, , ] <- t(Conj(temp)) %*% temp
        vect[f] <- fj * (index[j + 1] - index[j])
    }
    g <- colSums(sum_xwav, 1)/n

    r <- sum(log(eigen(g)$values)) - 2 * sum(d) * log(2) * sum(vect)/n
    return(r)
}
