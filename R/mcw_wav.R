mcw_wav <- function(xwav, index,  psih, grid_K, LU = NULL) 
{
    xwav <- as.matrix(xwav)
    k <- dim(xwav)[2]
    d_univ <- rep(0, k)
    for (ll in seq(1, k, 1)) {
        d_univ[ll] <- optimize(f = function(d){
            mcw_wav_eval(d, xwav = xwav[, ll], index = index, 
                LU = LU)
        }, lower = -10, upper = 10)$minimum
    }
    md <- d_univ
    if (k > 1) {
        md <- nlm(f = function(d){
            mcw_wav_eval(d, xwav = xwav, index = index, LU = LU)
        }, d_univ)$estimate
    }
    mg <- mcw_wav_cov_eval(md, xwav, index,  psih, grid_K, LU)
    list(d = md, cov = mg)
}

