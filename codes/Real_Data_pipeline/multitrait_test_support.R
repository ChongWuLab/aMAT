SumSqU <- function(U, CovS) {
    if (is.null(dim(CovS))) { # only one-dim:
        Tscore <- sum(U^2 / CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) Tscore <- 0
        pTg1 <- as.numeric(1 - pchisq(Tscore, 1))
    }
    else {
        # it's possible U=0 and Cov(U)=0:
        if (all(abs(U) < 1e-20)) {
            pTg1 <- 1
        } else {
            Tg1 <- t(U) %*% U
            ## distr of Tg1 is sum of cr Chisq_1:
            cr <- eigen(CovS, only.values = TRUE)$values
            ## approximate the distri by alpha Chisq_d + beta:
            alpha1 <- sum(cr * cr * cr) / sum(cr * cr)
            beta1 <- sum(cr) - (sum(cr * cr)^2) / (sum(cr * cr * cr))
            d1 <- (sum(cr * cr)^3) / (sum(cr * cr * cr)^2)
            alpha1 <- as.double(alpha1)
            beta1 <- as.double(beta1)
            d1 <- as.double(d1)
            pTg1 <- as.numeric(1 - pchisq((Tg1 - beta1) / alpha1, d1))
        }
    }
    return(pTg1)
}


SumSqU_fast <- function(U, beta1, alpha1, d1) {
    Tg1 <- t(U) %*% U

    pTg1 <- as.numeric(1 - pchisq((Tg1 - beta1) / alpha1, d1))
    return(pTg1)
}


MultiXcan <- function(Z.score, trait.cor, cutoff = 30) {
    Z.score <- as.numeric(Z.score)

    tmp <- svd(trait.cor)
    eigenvalue <- tmp$d
    eigenvalue[eigenvalue < max(eigenvalue) / cutoff] <- 0
    eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]

    multi.df <- sum(eigenvalue != 0)
    D <- diag(eigenvalue)
    t.mat2 <- tmp$u %*% D %*% t(tmp$v)

    multi.stat <- Z.score %*% t.mat2 %*% Z.score
    multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)


    return(list(multi.stat = multi.stat, multi.p = multi.p, multi.df = multi.df))
}


MultiXcan_fast <- function(Z.score, t.mat2, multi.df) { # tmp is  svd(trait.cor)
    multi.stat <- Z.score %*% t.mat2 %*% Z.score
    multi.p <- pchisq(multi.stat, df = multi.df, lower.tail = F)

    return(list(multi.stat = multi.stat, multi.p = multi.p, multi.df = multi.df))
}


########## SumTest########################
Sum <- function(U, CovS) {
    # it's possible U=0 and Cov(U)=0:
    if (all(abs(sum(U)) < 1e-20)) {
        pTsum <- 1
    } else {
        a <- rep(1, length(U))
        Tsum <- sum(U) / (sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric(1 - pchisq(Tsum^2, 1))
    }
    pTsum
}

########## SumTest########################
Sum_fast <- function(U, denom) {
    # it's possible U=0 and Cov(U)=0:

    Tsum <- sum(U) / denom
    pTsum <- as.numeric(1 - pchisq(Tsum^2, 1))

    pTsum
}
########## UminP Test########################
UminPd <- function(U, CovS) {
    if (is.null(dim(CovS))) { # only one-dim:
        Tu <- sum(U^2 / CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) Tu <- 0
        pTu <- as.numeric(1 - pchisq(Tu, 1))
    }
    else {
        pTu <- as.numeric(PowerUniv(U, CovS))
    }

    pTu
}


PowerUniv <- function(U, V) {
    n <- dim(V)[1]

    x <- as.numeric(max(abs(U)))
    TER <- as.numeric(1 - pmvnorm(lower = c(rep(-x, n)), upper = c(rep(x, n)), mean = c(rep(0, n)), sigma = V))

    return(TER)
}

# revise
omnibus_cor <- function(CovS, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.matall, multi.dfall, n.perm = 10000) {
    eV <- eigen(trait.cor)
    eV$values[eV$values < 0] <- 0

    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))

    T1s <- matrix(NA, n.perm, 6)
    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        pSum <- Sum_fast(tmp.z, sum_denom)

        pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
        #
        pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p

        T1s[i, ] <- c(pSum, pSSU, pAT1, pAT10, pAT30, pATall)
    }
    T1s[T1s > 0.99] <- 0.99
    T1s <- qnorm(1 - T1s)
    setbased_corEst <- cor(T1s)

    return(setbased_corEst)
}


omnibus_cor2 <- function(CovS, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall, n.perm = 10000) {
    eV <- eigen(trait.cor)
    eV$values[eV$values < 0] <- 0

    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))

    T1s <- matrix(NA, n.perm, 7)
    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        pSum <- Sum_fast(tmp.z, sum_denom)

        pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
        #
        pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pAT50 <- MultiXcan_fast(tmp.z, t.mat50, multi.df50)$multi.p

        pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p

        T1s[i, ] <- c(pSum, pSSU, pAT1, pAT10, pAT30, pAT50, pATall)
    }
    T1s[T1s > 0.99] <- 0.99
    T1s <- qnorm(1 - T1s)
    setbased_corEst <- cor(T1s)

    return(setbased_corEst)
}

# remove the extremely small singular values.
omnibus_cor3 <- function(CovS, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall, n.perm = 10000) {
    tmp <- svd(trait.cor)
    eigenvalue <- tmp$d
    #eigenvalue[eigenvalue < max(eigenvalue) / 10000] <- 0
    eigenvalue[eigenvalue < 0.001] <- 0
    D <- diag(eigenvalue)
    CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)

    T1s <- matrix(NA, n.perm, 7)
    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        pSum <- Sum_fast(tmp.z, sum_denom)

        pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
        #
        pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pAT50 <- MultiXcan_fast(tmp.z, t.mat50, multi.df50)$multi.p

        pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p

        T1s[i, ] <- c(pSum, pSSU, pAT1, pAT10, pAT30, pAT50, pATall)
    }
    T1s[T1s > 0.99] <- 0.99
    T1s <- qnorm(1 - T1s)
    setbased_corEst <- cor(T1s)

    return(setbased_corEst)
}


omnibus_with_minp_cor <- function(CovS, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.matall, multi.dfall, n.perm = 1000) {
    eV <- eigen(trait.cor)
    eV$values[eV$values < 0] <- 0

    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))

    T1s <- matrix(NA, n.perm, 7)
    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        pSum <- Sum_fast(tmp.z, sum_denom)

        pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
        #
        pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p
        minP <- PowerUniv(tmp.z, CovS)

        T1s[i, ] <- c(pSum, pSSU, pAT1, pAT10, pAT30, pATall, minP)
    }
    T1s[T1s > 0.99] <- 0.99
    T1s <- qnorm(1 - T1s)
    setbased_corEst <- cor(T1s)

    return(setbased_corEst)
}

