mcw_wav_cov_eval <-
function (d, xwav, index, psih, grid_K, LU) 
{
    if (is.matrix(xwav)) {
        N <- dim(xwav)[1]
        k <- dim(xwav)[2]
    }
    else {
        N <- length(xwav)
        k <- 1
    }
    xwav <- as.matrix(xwav, N, k)
    Jmax <- length(index) - 1
    if (is.null(LU)) {
        LU <- c(1, Jmax)
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
            temp <- xwav[(index[j] + 1):index[j + 1], ] %*% diag(2^(fj*d))
        }
        else {
            temp <- xwav[(index[j] + 1):index[j + 1]] * 2^(fj*d)
        }
        sum_xwav[f, , ] <- t(Conj(temp)) %*% (temp)
    }
    g_temp <- colSums(sum_xwav, 1)/n
    g <- matrix(0, k, k)
    g <- g_temp

    K <- K_eval(psih, grid_K, d)
    K <- Re(K)

    g <- g * (K^-1)

    return(g)
}

