K_eval <- function(psi_hat, u, d)  
{

    l <- length(psi_hat)
    k <- length(d)
    K_fct <- matrix(0, k, l)
    for (j in 1:k) {
        K_fct[j, ] <- (abs(u + (u == 0))^(-d[j]) - (u == 0)) * 
            psi_hat
    }
    K <- K_fct %*% t(Conj(K_fct)) * (max(u) - min(u))/l
    K <- Re(K)

    return(K)
}

