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


MAT <- function(Z.score, trait.cor, cutoff = 30) {
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


MAT_fast <- function(Z.score, t.mat2, multi.df) { # tmp is  svd(trait.cor)
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

aMAT <- function(Z, trait.cor,n.perm = 10000) {

    if(!(is.matrix(Z) | is.data.frame(Z))) {
        stop("Z (GWAS summary results) has to be a matrix or a data frame" )
    }
    if(sum(is.na(Z))!=0 ) {
        warning("Some missing values appear in Z. Impute with zero. Please double check.")
        Z[is.na(Z)] = 0
    }
    
    Z = apply(Z,2,as.numeric)
    
    n.snp = dim(Z)[1]
    res = matrix(NA,n.snp,5)
    rownames(res) = rownames(Z)
    colnames(res) = c("MAT(1)","MAT(10)","MAT(30)","MAT(50)","aMAT")
    
    tmp <- svd(trait.cor)
    eigenvalue <- tmp$d
    
    tmp <- svd(trait.cor)

    cat("Preparing intermidate statistics\n")

    eigenvalue <- tmp$d
    eigenvalue[eigenvalue < max(eigenvalue) / 1] <- 0
    eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
    multi.df1 <- sum(eigenvalue != 0)
    D <- diag(eigenvalue)
    t.mat1 <- tmp$u %*% D %*% t(tmp$v)
    invhalf.mat1 <- tmp$u %*% sqrt(D) %*% t(tmp$v)

    eigenvalue <- tmp$d
    eigenvalue[eigenvalue < max(eigenvalue) / 10] <- 0
    eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
    multi.df10 <- sum(eigenvalue != 0)
    D <- diag(eigenvalue)
    t.mat10 <- tmp$u %*% D %*% t(tmp$v)
    invhalf.mat10 <- tmp$u %*% sqrt(D) %*% t(tmp$v)


    eigenvalue <- tmp$d
    eigenvalue[eigenvalue < max(eigenvalue) / 30] <- 0
    eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
    multi.df30 <- sum(eigenvalue != 0)
    D <- diag(eigenvalue)
    t.mat30 <- tmp$u %*% D %*% t(tmp$v)
    invhalf.mat30 <- tmp$u %*% sqrt(D) %*% t(tmp$v)


    eigenvalue <- tmp$d
    eigenvalue[eigenvalue < max(eigenvalue) / 50] <- 0
    eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
    multi.df50 <- sum(eigenvalue != 0)
    D <- diag(eigenvalue)
    t.mat50 <- tmp$u %*% D %*% t(tmp$v)
    invhalf.mat50 <- tmp$u %*% sqrt(D) %*% t(tmp$v)


    tmp <- svd(trait.cor)
    eigenvalue <- tmp$d
    D <- diag(eigenvalue)
    CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)
    
    T1s <- matrix(NA, n.perm, 4)
    start.time <- proc.time()[3]

    cat("Estimating the correlation among MATs\n")

    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        pAT1 <- MAT_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MAT_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MAT_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pAT50 <- MAT_fast(tmp.z, t.mat50, multi.df50)$multi.p

        T1s[i, ] <- c(pAT1, pAT10, pAT30, pAT50)
    }

    colSums(T1s < 0.01,na.rm=T) / n.perm
    T1s.save <- T1s
    T1s[T1s > 0.99] <- 0.99
    T1s <- qnorm(1 - T1s)
    setbased_corEst1 <- cor(T1s)


    
    cat("Starting calculating p values\n")
    for (i in 1:n.snp) {
        tmp.z <- Z[i,]
        
        pAT1 <- MAT_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MAT_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MAT_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pAT50 <- MAT_fast(tmp.z, t.mat50, multi.df50)$multi.p

        omni_stat <- min(pAT1, pAT10, pAT30, pAT50)
        aMAT<- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(qnorm(1 - omni_stat), 4), sigma = setbased_corEst1)[1]

        res[i, ] <- c(pAT1, pAT10, pAT30, pAT50, aMAT)
    }
    
    return(res)
}
