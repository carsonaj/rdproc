# functions for collapsed gibbs sampler for sampling S | Y:

c_cgs <- function(y, al, a_mu, kap, a_tau, b_tau, burnin, num_samples) {

    # initialize sampler
    n <- length(y)
    s <- rep(1,n)
    samples <- vector("list", num_samples)

    # burnin
    for (t in 1:burnin) {
        s <- c_update_s(s, y, al, a_mu, kap, a_tau, b_tau)
    }

    # sample
    for (t in 1:num_samples) {
        s <- c_update_s(s, y, al, a_mu, kap, a_tau, b_tau)
        samples[[i]] <- s
        print("sample s is")
        print(s)
    }

    return(samples)
}

c_update_s <- function(s, y, al, a_mu, kap, a_tau, b_tau) {
    n <- length(s)
    c_check_length(y, n)
    samp <- s
    for (i in 1:n) {
        p <- c_get_full_s_i_prob(y, samp, i, al, a_mu, kap, a_tau, b_tau)
        print("p is")
        print(p)
        k <- length(p) - 1
        print("k is")
        print(k)
        upd_s_i <- sample(1:(k+1), 1, replace = T, prob = p)
        samp[i] <- upd_s_i
        print("samp is")
        print(samp)
    }

    return(samp)
}

# full conditional s_i | s_ni, y
c_get_full_s_i_prob <- function(y, s, i, al, a_mu, kap, a_tau, b_tau) {
    n <- length(y)
    c_check_length(s, n)
    ind_ni <- which(1:n != i)
    s_ni <- s[ind_ni]
    y_i <- y[i]
    y_ni <- y[ind_ni]
    k <- length(unique(s_ni))
    part_s_i_prob <- c_get_part_s_i_prob(n, s_ni, al)
    part_y_i_prob <- c_get_part_y_i_prob(y_i, y_ni, s_ni, a_mu, kap, a_tau, b_tau)

    c_check_length(part_y_i_prob, k+1)
    c_check_length(part_s_i_prob, k+1)

    full_si_prob <- part_s_i_prob * part_y_i_prob
    full_si_prob <- full_si_prob / sum(full_si_prob)

    c_check_length(full_si_prob, k+1)

    return(full_si_prob)
}

# partial conditional s_i | s_ni
c_get_part_s_i_prob <- function(n, s_ni, al) {
    k <- length(unique(s_ni))
    n_ni <- rep(0, k)
    s_i_prob <- rep(0, k+1)
    for (j in 1:k) {
        n_ni[j] <- length(which(s_ni == j))
        s_i_prob[j] <- n_ni[j] / (al + n - 1)
    }
    s_i_prob[k+1] <- al / (al + n - 1)
    s_i_prob <- s_i_prob / sum(s_i_prob)

    return(s_i_prob)
}

# partial conditional y_i|s_i=j, y_ni_j
c_get_part_y_i_prob <- function(y_i, y_ni, s_ni, a_mu, kap, a_tau, b_tau) {
    k <- length(unique(s_ni))
    y_i_prob <- rep(0, k+1)

    for (j in 1:k) {
        y_ni_j <- y_ni[which(s_ni == j)]
        y_ldata <- c(y_ni_j, y_i)
        y_i_prob[j] <- c_get_marg_y(y_ldata, kap, a_mu, a_tau, b_tau)
    }
    y_i_prob[k+1] <- c_get_marg_y(y_i, kap, a_mu, a_tau, b_tau)
    y_i_prob <- y_i_prob / sum(y_i_prob)

    return(y_i_prob)
}

c_check_length <- function(s, n) {
    if (length(s) != n) {
        stop("length of first arg must equal second arg")
    }
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# functions for sampling S_i | S_ni:

# return P(s_i = j|s_ni) unnormalized
c_get_p_si <- function(s, i, al) {
    n <- length(s)
    ind <- 1:n
    ind_ni <- ind[ind != i]
    s_ni <- s[ind_ni]
    k <- length(unique(s_ni))
    p_unnorm <- rep(0, k+1)
    n_ni <- rep(0, k)
    for (j in 1:k) {
        n_ni[j] <- length(s_ni == j)
        p_unnorm[j] <- n_ni[j] / (al + n -1)
    }
    p_unnorm[k+1] <- al / (al + n - 1)

    return(p_unnorm)
}
#---------------------------------------------------------------------------------------

# return P(s_{n+1} = j|s) normalized
c_get_p_snp1 <- function(s, al) {
    n <- length(s)
    k <- length(unique(s))
    p_norm <- rep(0, k+1)
    n_nnp1 <- rep(0, k)
    for (j in 1:k) {
        n_nnp1[j] <- length(s == j)
        p_norm[j] <- n_nnp1[j] / (al + n -1)
    }
    p_norm[k+1] <- al / (al + n - 1)

    #normalize
    p_norm <- p / sum(p)

    return(p_norm)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# functions for norm-norm-gamma posterior
# y_i | mu, tau ~ N(mu, tau^(-1))
# mu | tau ~ N(a_mu, k * tau^(-1))
# tau ~ gamma(a_tau, b_tau)

c_get_mu_n <- function(n, kap_inv, ybar, a_mu) {
    mu_n <- (kap_inv * a_mu + n * ybar) / (n + kap_inv)

    return(mu_n)
}
#---------------------------------------------------------------------------------------

c_get_al_n <- function(n, a_tau) {
    al_n <- n / 2 + a_tau

    return(al_n)
}
#---------------------------------------------------------------------------------------

c_get_bet_n <- function(n, kap_inv, ybar, s_sq, a_mu, b_tau, mu_n) {
    bet_n <- .5 * (n * (s_sq + ybar^2) + kap_inv * a_mu^2 - (n + kap_inv) * mu_n^2)

    return(bet_n)
}
#---------------------------------------------------------------------------------------

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
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------






