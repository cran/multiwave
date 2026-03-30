mcw <- function(x,filter_complex,LU=NULL,J=10){

  ## Computes the multivariate wavelet Whittle estimation for the of
  ## the long-memory parameter d and the long-run covariance matrix.
  ##
  ## 	INPUT
  ##    x	 Data (n*k vector)
  ##		filter_complex	Complex or real wavelet filter
  ##    LU	bivariate vector (optional) containing
  ##			  L, the lowest resolution in wavelet decomposition
  ##        U, the maximal resolution in wavelet decomposition
  ##    J  2^J corresponds to the size of the grid for the discretisation of the wavelet.
  ##       The default value is set to 10.
  ##
  ##	OUTPUT
  ##    d	 Long-range parameter estimation
  ##		cov	Long-run covariance matrix estimation
  ##
  ##						Achard & Gannaz (2024)
  ##______________________________________________________________________________



        N <- dim(x)[1]
        k <- dim(x)[2]

     xwav <- matrix(0, N, k)
    for (j in 1:k) {
        xx <- x[, j]
        if( is.list(filter_complex) ){
          resw <- DWTcomplex(xx, filter_complex)
          res_complex_psi <- psi_hat_exact_complex(filter_complex$h,filter_complex$g,J)
          psih_complex <- res_complex_psi$psih
          grid_complex <- res_complex_psi$grid
        }else{
          resw <- DWTexact(xx, filter_complex)
          res_psi <- psi_hat_exact(filter_complex,J)
          psih_complex <- res_psi$psih
          grid_complex <- res_psi$grid
        }

        xwav_temp <- resw$dwt
        index <- resw$indmaxband
        Jmax <- resw$Jmax
        xwav[1:index[Jmax], j] <- xwav_temp
    }
    new_xwav <- matrix(0, min(index[Jmax], N), k)
    if (index[Jmax] < N) {
        new_xwav[(1:(index[Jmax])), ] <- xwav[(1:(index[Jmax])), ]
    }
    xwav <- new_xwav
    index <- c(0, index)

    res <- mcw_wav(xwav,index,psih_complex,grid_complex,LU)

    return(res)


}

