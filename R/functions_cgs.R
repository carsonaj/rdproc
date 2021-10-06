# return P(s_i = j|s_ni)
c_cond_den_s_i_given_s_ni <- function(s, i, al) {
    n <- length(s)
    ind <- 1:n
    ind_ni <- ind[ind != i]
    s_ni <- s[ind_ni]
    k <- length(unique(s_ni))
    p <- rep(0, k+1)
    n_ni <- rep(0, k)
    for (j in 1:k) {
        n_ni[j] <- length(s_ni == j)
        p[j] <- n_ni[j] / (al + n -1)
    }
    p[k+1] <- al / (al + n - 1)

    return(p)
}

# y_i | mu, tau ~ N(mu, tau^(-1))
# mu | tau ~ N(a_mu, k * tau^(-1))
# tau ~ gamma(a_tau, b_tau)
c_get_mu_n <- function(n, kap_inv, ybar, a_mu) {
    mu_n <- (kap_inv * a_mu + n * ybar) / (n + kap_inv)

    return(mu_n)
}

c_get_al_n <- function(n, a_tau) {
    al_n <- n / 2 + a_tau

    return(al_n)
}

c_get_bet_n <- function(n, kap_inv, ybar, s_sq, a_mu, b_tau, mu_n) {
    bet_n <- .5 * (n * (s_sq + ybar^2) + kap_inv * a_mu - (n + kap_inv) * mu_n^2)

    return(bet_n)
}

c_get_marg_y <- function(y, kap, a_mu, a_tau, b_tau) {
    n <- length(y)
    kap_inv <- kap^(-1)
    ybar <- mean(y)
    s_sq <- mean((y - ybar)^2)
    mu_n <- c_get_mu_n(n, kap_inv, ybar, a_mu)
    al_n <- c_get_al_n(n, a_tau)
    bet_n <- c_get_bet_n(n, kap_inv, ybar, s_sq, a_mu, b_tau, mu_n)

    marg_y <- (2 * pi)^(- n / 2) * (1 + kap * n)^(1 / 2) *
        b_tau * gamma(al_n) / (bet_n^(al_n) * gamma(a_tau))

    return(marg_y)
}






