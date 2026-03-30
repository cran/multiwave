psi_hat_exact <- function(filter, J = 10)
{

  ## Computes the discrete Fourier transform of the real wavelet associated to
  ## the given filter
  ##
  ## The wavelet is given by the function makescalingfunction of Fay et al
  ## (2009).
  ##
  ## The length of the Fourier transform is equal to the length of the grid
  ## where the wavelet is evaluated
  ##
  ## 	INPUT	filter  Quadratic mirror filter of the wavelet
  ##              J	Number of scales where the wavelet is evaluated
  ##
  ##	OUTPUT	fct	Values of the discrete Fourier transform of the wavelet
  ##		grid	Frequencies where the Fourier transform is evaluated
  ##
  ##                                           		Achard & Gannaz (2014)
  ##________________________________________________________________________________

    Jmax <- J - 3

    res <- scaling_function(filter, J)

    phi <- res$phi
    psi <- res$psi
    Ll <- length(psi)
    psi <- psi/sqrt(2^J)

    Nfft <- 2^ceiling(log(Ll)/log(2))
    Te <- (length(filter) - 1) * 2^Jmax
    psi <- c(psi, rep(0, (Nfft - Ll)))
    psih <- fft(psi,inverse=TRUE)/sqrt(pi * Te)
    grid <- pi * Te * c(seq((1/Nfft), 1/2, (1/Nfft)), seq((-1/2),
        (-1/Nfft), (1/Nfft)))
	new_grid<-sort(grid,index.return=T)
	final_grid<-new_grid$x
	final_psih<-psih[new_grid$ix]

    phih <- fft(phi)/sqrt(pi * Te)
	final_phih<-phih[new_grid$ix]

    final_psih <- final_psih
    final_grid <- final_grid

    list(phih = final_phih, psih = final_psih, grid = final_grid)
}



psi_hat_exact_complex <- function(h, g, J = 10)#, complex=F)
{

  ## Computes the discrete Fourier transform of the common-factor complex
  ## wavelet associated to the given filters (h,g) obtained with the hwlet function
  ##
  ## The wavelet is given by the function makescalingfunction of Fay et al
  ## (2009).
  ##
  ## The length of the Fourier transform is equal to the length of the grid
  ## where the wavelet is evaluated
  ##
  ## 	INPUT	filter  Filters of a complex Common-Factor wavelet
  ##              J	Number of scales where the wavelet is evaluated
  ##
  ##	OUTPUT	fct	Values of the discrete Fourier transform of the wavelet
  ##		grid	Frequencies where the Fourier transform is evaluated
  ##
  ##                                           		Achard & Gannaz (2024)
  ##________________________________________________________________________________

    Jmax <- J - 3

    #####################
    res <- scaling_function(h, J)

    phi <- res$phi
    psi <- res$psi
    Ll <- length(psi)
    psi <- psi/sqrt(2^J)

    Nfft <- 2^ceiling(log(Ll)/log(2))
    Te <- (length(h) - 1) * 2^Jmax
    psi <- c(psi, rep(0, (Nfft - Ll)))
    psih <- fft(psi,inverse=TRUE)/sqrt(pi * Te)
    grid <- pi * Te * c(seq((1/Nfft), 1/2, (1/Nfft)), seq((-1/2),
        (-1/Nfft), (1/Nfft)))
	new_grid<-sort(grid,index.return=T)
	final_grid<-new_grid$x
	final_psih<-psih[new_grid$ix]

    phi <- phi/sqrt(2^J)
    phi <- c(phi, rep(0, (Nfft - Ll)))
    phih <- fft(phi)/sqrt(pi * Te)
	final_phih<-phih[new_grid$ix]

    ######################
    res <- scaling_function(g, J)

    phi <- res$phi
    psi <- res$psi
    psi <- psi/sqrt(2^J)
    psi <- c(psi, rep(0, (Nfft - Ll)))
    psih <- fft(psi,inverse=TRUE)/sqrt(pi * Te)
	final_psih<-final_psih+1i*psih[new_grid$ix]

    phi <- phi/sqrt(2^J)
    phi <- c(phi, rep(0, (Nfft - Ll)))
    phih <- fft(phi)/sqrt(pi * Te)
	final_phih <- final_phih+1i*phih[new_grid$ix]

    ######################
    final_phih <- final_phih/sqrt(2)
    final_psih <- final_psih/sqrt(2)
    final_grid <- final_grid

    list(phih = final_phih, psih = final_psih, grid = final_grid)
}
