test_that("functions_bgs.R works", {
    source("../../R/functions_bgs.R")
    #source("R/functions_bgs.R")

    # tests:

    # test for bgs:
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    k <- 10
    al <- 1
    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1
    burnin <- 10
    num_samples <- 10
    samples <- bgs(y, k, al, a_tau, b_tau, a_mu, kap, burnin, num_samples)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_check_length:
    z1 <- 1:5
    y1 <- seq(1, 2, length.out = 5)
    z2 <- 1:6
    expect_error(b_check_length(z1, y1), NA)
    expect_error(b_check_length(z2, y1),
                 "length of first arg must equal length of second arg",
                 fixed = T)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_check_normalized:
    p <- c(5.2, 1.333, 3.8393984839)
    expect_error(b_check_normalized(p), "argument must be normalized", fixed = T)

    p <- p / sum(p)
    expect_error(b_check_normalized(p), NA)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_z_info:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- b_get_z_info(z, y)
    expect_equal(z_info$z_vals, c(2, 5, 7, 9))
    expect_equal(z_info$grps[[1]], list("y_j" = c(1.2, -3.5), "n_j" = 2, "ybar_j" = -1.15))
    expect_equal(z_info$grps[[4]], list("y_j" = -1, "n_j" = 1, "ybar_j" = -1))
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_n:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- b_get_z_info(z, y)

    n <- b_get_n(10, z_info)
    n_check <- c(0, 2, 0, 0, 2, 0, 1, 0, 1, 0)

    expect_equal(n, n_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_tau_params:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- b_get_z_info(z, y)

    tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    a1 <- tau_params[[1]]$a
    b1 <- tau_params[[1]]$b
    a4 <- tau_params[[4]]$a
    b4 <- tau_params[[4]]$b

    kap_inv <- kap^(1)
    n_1 <- z_info$grps[[1]]$n_j
    y_1 <- z_info$grps[[1]]$y_j
    ybar_1 <- z_info$grps[[1]]$ybar_j
    n_4 <- z_info$grps[[4]]$n_j
    y_4 <- z_info$grps[[4]]$y_j
    ybar_4 <- z_info$grps[[4]]$ybar_j


    a1_check <- a_tau + n_1 / 2
    b1_check <- (sum((y_1 - ybar_1)^2) +
                     (kap_inv * n_1 * (ybar_1 - a_mu)^2) / (2 * (kap_inv + n_1))) / 2
    a4_check <- a_tau + n_4 / 2
    b4_check <- (sum((y_4 - ybar_4)^2) +
                     (kap_inv * n_4 * (ybar_4 - a_mu)^2) / (2 * (kap_inv + n_4))) / 2

    expect_equal(a1, a1_check)
    expect_equal(b1, b1_check)
    expect_equal(a4, a4_check)
    expect_equal(b4, b4_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_tau:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- b_get_z_info(z, y)
    tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- b_get_tau(k, z_info, tau_params, a_tau, b_tau)

    tau_1 <- tau[1]
    tau_2 <- tau[2]

    set.seed(1000)
    a <- tau_params[[1]]$a
    b <- tau_params[[1]]$b
    tau_1_check <- rgamma(1, a_tau, b_tau)
    tau_2_check <- rgamma(1, a, b)

    expect_equal(tau_1, tau_1_check)
    expect_equal(tau_2, tau_2_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_mu_params:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- b_get_z_info(z, y)
    tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- b_get_tau(k, z_info, tau_params, a_tau, b_tau)
    mu_params <- b_get_mu_params(z_info, a_mu, kap)

    set.seed(1000)
    a1 <- mu_params[[1]]$a
    b1 <- mu_params[[1]]$b
    a4 <- mu_params[[4]]$a
    b4 <- mu_params[[4]]$b

    kap_inv <- kap^(-1)
    n_1 <- z_info$grps[[1]]$n_j
    y_1 <- z_info$grps[[1]]$y_j
    ybar_1 <- z_info$grps[[1]]$ybar_j
    n_4 <- z_info$grps[[4]]$n_j
    y_4 <- z_info$grps[[4]]$y_j
    ybar_4 <- z_info$grps[[4]]$ybar_j


    a1_check <- (kap_inv * a_mu + n_1 * ybar_1) / (kap_inv + n_1)
    b1_check <- (kap_inv + n_1)^(-1)
    a4_check <- (kap_inv * a_mu + n_4 * ybar_4) / (kap_inv + n_4)
    b4_check <- (kap_inv + n_4)^(-1)

    expect_equal(a1, a1_check)
    expect_equal(b1, b1_check)
    expect_equal(a4, a4_check)
    expect_equal(b4, b4_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_mu:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- b_get_z_info(z, y)
    tau_params <- b_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- b_get_tau(k, z_info, tau_params, a_tau, b_tau)
    mu_params <- b_get_mu_params(z_info, a_mu, kap)

    expect_error(b_get_mu(k, z_info, mu_params, kap, tau, a_mu), NA)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------


    # test for b_get_pi:
    nu <- c(.25, .25, .5)
    pi <- b_get_pi(nu)
    pi_check <- c(.25, .25 * .75, .5 * .75^2, .5 * .75^2)

    expect_equal(pi, pi_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_z
    y <- c(-.25, 0, .25)
    pi <- c(.3, .3, .4)
    mu <- c(-1, 0, 1)
    tau <- c(1, 1, 1)

    z <- b_get_z(y, pi, mu, tau)

    expect_equal(length(z), 3)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for b_get_nu:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- b_get_z_info(z, y)
    k <- 10
    al <- 1

    nu <- b_get_nu(k, al, z_info)

    expect_equal(length(nu), k-1)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
})

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

test_that("functions_fmmgs.R works", {
    source("../../R/functions_fmmgs.R")
    #source("R/functions_fmmgs.R")

    # tests:

    # test for fmmgs:
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    k <- 10
    al <- 1
    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1
    burnin <- 10
    num_samples <- 10
    samples <- fmmgs(y, k, al, a_tau, b_tau, a_mu, kap, burnin, num_samples)

    # test for f_check_length:
    z1 <- 1:5
    y1 <- seq(1, 2, length.out = 5)
    z2 <- 1:6
    expect_error(f_check_length(z1, y1), NA)
    expect_error(f_check_length(z2, y1),
                 "length of first arg must equal length of second arg",
                 fixed = T)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_check_normalized:
    p <- c(5.2, 1.333, 3.8393984839)
    expect_error(f_check_normalized(p), "argument must be normalized", fixed = T)

    p <- p / sum(p)
    expect_error(f_check_normalized(p), NA)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_z_info:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- f_get_z_info(z, y)
    expect_equal(z_info$z_vals, c(2, 5, 7, 9))
    expect_equal(z_info$grps[[1]], list("y_j" = c(1.2, -3.5), "n_j" = 2, "ybar_j" = -1.15))
    expect_equal(z_info$grps[[4]], list("y_j" = -1, "n_j" = 1, "ybar_j" = -1))
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_n:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- f_get_z_info(z, y)

    n <- f_get_n(10, z_info)
    n_check <- c(0, 2, 0, 0, 2, 0, 1, 0, 1, 0)

    expect_equal(n, n_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_tau_params:
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- f_get_z_info(z, y)

    tau_params <- f_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    a1 <- tau_params[[1]]$a
    b1 <- tau_params[[1]]$b
    a4 <- tau_params[[4]]$a
    b4 <- tau_params[[4]]$b

    kap_inv <- kap^(1)
    n_1 <- z_info$grps[[1]]$n_j
    y_1 <- z_info$grps[[1]]$y_j
    ybar_1 <- z_info$grps[[1]]$ybar_j
    n_4 <- z_info$grps[[4]]$n_j
    y_4 <- z_info$grps[[4]]$y_j
    ybar_4 <- z_info$grps[[4]]$ybar_j


    a1_check <- a_tau + n_1 / 2
    b1_check <- (sum((y_1 - ybar_1)^2) +
                     (kap_inv * n_1 * (ybar_1 - a_mu)^2) / (2 * (kap_inv + n_1))) / 2
    a4_check <- a_tau + n_4 / 2
    b4_check <- (sum((y_4 - ybar_4)^2) +
                     (kap_inv * n_4 * (ybar_4 - a_mu)^2) / (2 * (kap_inv + n_4))) / 2

    expect_equal(a1, a1_check)
    expect_equal(b1, b1_check)
    expect_equal(a4, a4_check)
    expect_equal(b4, b4_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_tau:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- f_get_z_info(z, y)
    tau_params <- f_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- f_get_tau(k, z_info, tau_params, a_tau, b_tau)

    tau_1 <- tau[1]
    tau_2 <- tau[2]

    set.seed(1000)
    a <- tau_params[[1]]$a
    b <- tau_params[[1]]$b
    tau_1_check <- rgamma(1, a_tau, b_tau)
    tau_2_check <- rgamma(1, a, b)

    expect_equal(tau_1, tau_1_check)
    expect_equal(tau_2, tau_2_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_mu_params:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- f_get_z_info(z, y)
    tau_params <- f_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- f_get_tau(k, z_info, tau_params, a_tau, b_tau)
    mu_params <- f_get_mu_params(z_info, a_mu, kap)

    set.seed(1000)
    a1 <- mu_params[[1]]$a
    b1 <- mu_params[[1]]$b
    a4 <- mu_params[[4]]$a
    b4 <- mu_params[[4]]$b

    kap_inv <- kap^(-1)
    n_1 <- z_info$grps[[1]]$n_j
    y_1 <- z_info$grps[[1]]$y_j
    ybar_1 <- z_info$grps[[1]]$ybar_j
    n_4 <- z_info$grps[[4]]$n_j
    y_4 <- z_info$grps[[4]]$y_j
    ybar_4 <- z_info$grps[[4]]$ybar_j


    a1_check <- (kap_inv * a_mu + n_1 * ybar_1) / (kap_inv + n_1)
    b1_check <- (kap_inv + n_1)^(-1)
    a4_check <- (kap_inv * a_mu + n_4 * ybar_4) / (kap_inv + n_4)
    b4_check <- (kap_inv + n_4)^(-1)

    expect_equal(a1, a1_check)
    expect_equal(b1, b1_check)
    expect_equal(a4, a4_check)
    expect_equal(b4, b4_check)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_mu:
    set.seed(1000)
    k <- 10
    z <- c(2, 2, 5, 5, 7, 9)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)

    a_tau <- 1
    b_tau <- 1
    a_mu <- 0
    kap <- 1

    z_info <- f_get_z_info(z, y)
    tau_params <- f_get_tau_params(z_info, a_tau, b_tau, a_mu, kap)
    tau <- f_get_tau(k, z_info, tau_params, a_tau, b_tau)
    mu_params <- f_get_mu_params(z_info, a_mu, kap)

    expect_error(f_get_mu(k, z_info, mu_params, kap, tau, a_mu), NA)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------


    # test for f_get_pi:
    set.seed(1000)
    z <- c(1, 2, 2, 2, 2, 3)
    y <- c(1.2, -3.5, 0.4, 2, -.4, -1)
    z_info <- f_get_z_info(z, y)
    k <- 4
    al <- 1
    pi <- f_get_pi(k, z_info, al)

    expect_equal(length(pi), k)
    expect_equal(sum(pi[2] >= pi) == k, T)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------

    # test for f_get_z
    y <- c(-.25, 0, .25)
    pi <- c(.3, .3, .4)
    mu <- c(-1, 0, 1)
    tau <- c(1, 1, 1)

    z <- f_get_z(y, pi, mu, tau)

    expect_error(f_check_length(z, 1:3), NA)
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------
})



