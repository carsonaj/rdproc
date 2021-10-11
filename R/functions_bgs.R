# blocked gibbs sampler:

# y_i: data, i in [n]

# nu_l: l in [k-1]
# nu_l iid beta(1, al)

# pi_j: stick breaking, j in [k]
# for j in [k-1], pi_j = vu_j * prod (1-nu_l) where 1 <= l <= j-1
# for j = k, pi_j = prod (1-nu_l) where 1 <= l <= j-1

# z_i: latent variable for dealing with mixture, i in [n]
# z_i|nu in [k] are iid and P(z_i = j|nu) = pi_j

# tau_j: j in [k]
# tau_j iid gamma(a_tau, b_tau)

# mu_j: j in [k]
# mu_j|tau_j iid normal(a_mu, kap * tau_j^(-1))

# n_j: j in [k]
# n_j = #{i : z_i = j}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# posterior predictive density estimate
# yn is y_new, the new y value to be drawn given the data y
#' @export
bppd <- function(yn, y, samples) {
    k <- length(samples[[1]]$pi)
    N <- length(samples)

    yn_g_theta <- rep(0, N)
    for (i in 1:N) {
        yn_g_theta_i <- rep(0, k)
        for (j in 1:k) {
            pi <- samples[[i]]$pi
            mu <- samples[[i]]$mu
            tau <- samples[[i]]$tau
            yn_g_theta_i[j] = pi[j] * dnorm(yn, mu[j], tau[j]^(-.5))
        }

        yn_g_theta[i] <- sum(yn_g_theta_i)
    }

    yn_g_y <- mean(yn_g_theta)

    return(yn_g_y)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# blocked_gibbs_sampler: sample tau, mu, pi, z| y
#' @export
bgs <- function(y, k, al, a_tau, b_tau, a_mu, kap, burnin, num_samples) {
    samples <- vector("list", num_samples)
    # initialize parameters
    n <- length(y)
    z <- rep(1, n)
    z_info <- b_get_z_info(z, y)
    nu <- rep(.5, k-1)
    tau <- rep(1, k)
    mu <- rep(0, k)

    # burnin iterations
    for (i in 1:burnin) {
        nu <- b_get_nu(k, al, z_info)
        pi <- b_get_pi(nu)
        z <- b_get_z(y, pi, mu, tau)
        z_info <- b_get_z_info(z, y)
        tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
        tau <- b_get_tau(k, z_info, tau_params, a_tau, b_tau)
        mu_params <- b_get_mu_params(z_info, a_mu, kap)
        mu <- b_get_mu(k, z_info, mu_params, kap, tau, a_mu)
    }

    # sample iterations
    for (i in 1:num_samples) {
        nu <- b_get_nu(k, al, z_info)
        pi <- b_get_pi(nu)
        z <- b_get_z(y, pi, mu, tau)
        z_info <- b_get_z_info(z, y)
        tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
        tau <- b_get_tau(k, z_info, tau_params, a_tau, b_tau)
        mu_params <- b_get_mu_params(z_info, a_mu, kap)
        mu <- b_get_mu(k, z_info, mu_params, kap, tau, a_mu)

        samples[[i]] <- list("tau" = tau, "mu" = mu, "pi" = pi, "z" = z)
    }

    return(samples)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_check_length <- function(z, y) {
    if (length(z) != length(y)) {
        stop("length of first arg must equal length of second arg")
    }
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_check_normalized <- function(p) {
    if (abs(sum(p) - 1) > .001) {
        stop("argument must be normalized")
    }
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_z_info <- function(z, y) {
    b_check_length(z, y)
    z_vals <- unique(z)
    l <- length(z_vals)
    grps <- vector("list", l)

    for (j in 1:l) {
        ind <- which(z == z_vals[j])
        n_j <- length(ind)
        y_j <- y[ind]
        ybar_j <- mean(y_j)
        grps[[j]] <- list("y_j" = y_j, "n_j" = n_j, "ybar_j" = ybar_j)
    }

    info <- list("z_vals" = z_vals, "grps" = grps)

    return(info)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_n <- function(k, z_info) {
    n <- rep(0, k)
    z_vals <- z_info$z_vals
    grps <- z_info$grps

    for (j in 1:k) {
        if (j %in% z_vals){
            ind <- which(z_vals == j)
            n[j] <- grps[[ind]]$n_j
        } else {
            n[j] <- 0
        }
    }

    # if (sum(n) != length(y)) {
    #     stop("sum of n must equal length of y")
    # }

    return(n)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_tau_params <- function(z_info, a_tau, b_tau, a_mu, kap) {
    l <- length(z_info$z_vals)
    grps <- z_info$grps
    tau_params <- vector("list", l)
    kap_inv <- kap^(-1)
    for (j in 1:l) {
        n_j <- grps[[j]]$n_j
        y_j <- grps[[j]]$y_j
        ybar_j <- grps[[j]]$ybar_j

        a <- a_tau + n_j / 2
        b <- (sum((y_j - ybar_j)^2) +
                  (kap_inv * n_j * (ybar_j - a_mu)^2) / (2 * (kap_inv + n_j))) / 2

        tau_params[[j]] <- list("a" = a, "b" = b)
    }

    return(tau_params)
}

b_get_tau <- function(k, z_info, tau_params, a_tau, b_tau) {
    z_vals <- z_info$z_vals
    tau <- rep(0, k)
    mu <- rep(0, k)
    for (j in 1:k) {
        if (j %in% z_info$z_vals) {
            index <- which(z_vals == j)
            a <- tau_params[[index]]$a
            b <- tau_params[[index]]$b
            tau[j] <- rgamma(1, a, b)
        } else {
            tau[j] <- rgamma(1, a_tau, b_tau)
        }
    }

    return(tau)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_mu_params <- function(z_info, a_mu, kap) {
    l <- length(z_info$z_vals)
    grps <- z_info$grps
    mu_params <- vector("list", l)
    kap_inv <- kap^(-1)
    for (j in 1:l) {
        n_j <- grps[[j]]$n_j
        y_j <- grps[[j]]$y_j
        ybar_j <- grps[[j]]$ybar_j

        a <- (kap_inv * a_mu + n_j * ybar_j) / (kap_inv + n_j)
        b <- (kap_inv + n_j)^(-1)

        mu_params[[j]] <- list("a" = a, "b" = b)
    }

    return(mu_params)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_mu <- function(k, z_info, mu_params, kap, tau, a_mu) {
    z_vals <- z_info$z_vals
    mu <- rep(0, k)
    for (j in 1:k) {
        if (j %in% z_vals) {
            index <- which(z_vals == j)
            a <- mu_params[[index]]$a
            b <- mu_params[[index]]$b
            mu[j] <- rnorm(1, a, sqrt(b * tau[j]^(-1)))
        } else {
            mu[j] <- rnorm(1, a_mu, sqrt(kap * tau[j]^(-1)))
        }
    }

    return(mu)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_pi <- function(nu) {
    k <- length(nu) + 1
    pi <- rep(0, k)
    pi[1] <- nu[1]
    for (j in 2:(k-1)) {
        pi[j] <- nu[j] * prod(1 - nu[1:(j-1)])
    }
    pi[k] <- prod(1 - nu[1:(k-1)])

    return(pi)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_z <- function(y, pi, mu, tau) {
    b_check_length(pi, mu)
    b_check_length(pi, tau)
    k <- length(pi)
    n <- length(y)
    z <- rep(0, n)

    for (i in 1:n) {
        p <- rep(0, k)
        for (j in 1:k) {
            p = pi * dnorm(y[i], mu, sqrt(tau^(-1)))
        }
        #p[which(is.na(p))] <- 0
        p[which(is.infinite(p))] <- 0
        p <- p / sum(p)
        #f_check_normalized(p)
        p[which(is.na(p))] <- 0
        #print(p)
        z[i] <- sample(1:k, 1, replace = TRUE, prob = p)
    }

    return(z)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

b_get_nu <- function(k, al, z_info) {
    nu <- rep(0, k-1)
    n <- b_get_n(k, z_info)
    z_vals <- z_info$z_vals
    grps <- z_info$grps

    for (j in 1:(k-1)) {
        nu[j] <- rbeta(1, n[j] + 1, al + sum(n[(j+1):k]))
    }

    return(nu)
}
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
