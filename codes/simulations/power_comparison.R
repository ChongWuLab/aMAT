setwd("/gpfs/research/chongwu/Chong/UKBiobankImage")
library("MASS")
library("Matrix")
require(compiler)
enableJIT(4)

library(data.table)
library(mvtnorm)
library(MTAR)
source("/gpfs/research/chongwu/Chong/Multitrait-snp/multitrait_test_support.R")
source("/gpfs/research/chongwu/Chong/Multitrait-snp/FunctionSet.R")

suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))

input.cutoff = job
# pos <- readRDS("processed_hap_pos.rds")

id <- readRDS("IDP_id.rds")
# need to use 504 cores; for others, using 264 cores.
#idpgroup <- "fiveIDPs"
#setting <- "setting0"

Freesurfer.set <- c("Freesurfer", "Freesurfer_volume", "Freesurfer_area", "Freesurfer_thickness","fiveIDPs","twentyfiveIDPs")

if (idpgroup %in% Freesurfer.set) {
    if (idpgroup == "Freesurfer") {
        id <- id[grepl("volume_", id[, 2]) | grepl("_thickness", id[, 2]) | grepl("_area", id[, 2]), ]
        id <- id[8:dim(id)[1], ]
    }

    if (idpgroup == "Freesurfer_volume" | idpgroup== "fiveIDPs" | idpgroup== "twentyfiveIDPs" ) {
        id <- id[grepl("volume_", id[, 2]), ]
        id <- id[8:dim(id)[1], ]
    }

    if (idpgroup == "Freesurfer_thickness") {
        id <- id[grepl("_thickness", id[, 2]), ]
    }

    if (idpgroup == "Freesurfer_area") {
        id <- id[grepl("_area", id[, 2]), ]
    }

    remain <- readRDS("/gpfs/research/chongwu/Chong/Multitrait-snp/data/remain_Freesurfer_cor.rds")
    remain.tmp <- c(remain[, 1], remain[, 2])
    remain.tmp <- table(remain.tmp)
    remain.tmp <- sort(remain.tmp)
    remove.idp <- names(remain.tmp)[remain.tmp > 470]
    remove.idp <- gsub(".sumstats", "", remove.idp)

    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    id <- id[!id$usedIDP %in% remove.idp, ]

    trait.cor <- readRDS(paste("/gpfs/research/chongwu/Chong/Multitrait-snp/data/", "Freesurfer", "_cor.rds", sep = ""))

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}

trait.diag <- readRDS("/gpfs/research/chongwu/Chong/Multitrait-snp/ldsc/results/res_genetic_cor_diagonal.rds")
trait.diag[, 3] <- gsub(".sumstats", "", trait.diag[, 3])
rownames(trait.diag) <- trait.diag[, 3]

trait.diag <- trait.diag[rownames(trait.cor), ]
trait.diag[is.na(trait.diag[, 4]), 4] <- 1
diag(trait.cor) <- trait.diag[, 4]

if(idpgroup== "fiveIDPs") {
    trait.cor = trait.cor[1:5,1:5]
}

if(idpgroup== "twentyfiveIDPs") {
    trait.cor = trait.cor[1:25,1:25]
}

tmp <- svd(trait.cor)
eigenvalue <- tmp$d

D <- diag(eigenvalue)
trait.cor <- tmp$u %*% D %*% t(tmp$v)

# Run tests;
a <- rep(1, dim(trait.cor)[1])
sum_denom <- (sqrt(as.numeric(t(a) %*% trait.cor %*% (a))))

# SSU
cr <- eigen(trait.cor, only.values = TRUE)$values
## approximate the distri by alpha Chisq_d + beta:
alpha1 <- sum(cr * cr * cr) / sum(cr * cr)
beta1 <- sum(cr) - (sum(cr * cr)^2) / (sum(cr * cr * cr))
d1 <- (sum(cr * cr)^3) / (sum(cr * cr * cr)^2)
alpha1 <- as.double(alpha1)
beta1 <- as.double(beta1)
d1 <- as.double(d1)


tmp <- svd(trait.cor)

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


eigenvalue <- tmp$d
eigenvalue[eigenvalue < 0.001] <- 0
eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
multi.dfallrev <- sum(eigenvalue != 0)
D <- diag(eigenvalue)
t.matallrev <- tmp$u %*% D %*% t(tmp$v)
invhalf.matall <- tmp$u %*% sqrt(D) %*% t(tmp$v)


eigenvalue <- tmp$d
eigenvalue[eigenvalue < 0] <- 0
eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
multi.dfall <- sum(eigenvalue != 0)
D <- diag(eigenvalue)
t.matall <- ginv(trait.cor)
# t.matall <- tmp$u %*% D %*% t(tmp$v)
invhalf.matall <- tmp$u %*% sqrt(D) %*% t(tmp$v)


SampleSize <- rep(8411, dim(trait.cor)[1])

Wi <- matrix(SampleSize, nrow = 1)
sumW <- sqrt(sum(Wi^2))
W <- Wi / sumW
Sigma <- ginv(trait.cor)
# methods.cor <- omnibus_cor2(  trait.cor, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall,n.perm = 10000)


tmp <- svd(trait.cor)
eigenvalue <- tmp$d
# eigenvalue[eigenvalue < max(eigenvalue) / 10000] <- 0
eigenvalue[eigenvalue < 0] <- 0
D <- diag(eigenvalue)
CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)


cor.nperm <- 10000

T1s <- matrix(NA, cor.nperm, 4)
start.time <- proc.time()[3]

for (i in 1:cor.nperm) {
    tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
    tmp.z <- as.vector(tmp.z)

    #  pSum <- Sum_fast(tmp.z, sum_denom)
    # pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
    #
    pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
    pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
    pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
    pAT50 <- MultiXcan_fast(tmp.z, t.mat50, multi.df50)$multi.p

    T1s[i, ] <- c(pAT1, pAT10, pAT30, pAT50)
}

colSums(T1s < 0.01) / cor.nperm
T1s.save <- T1s
T1s[T1s > 0.99] <- 0.99
T1s <- qnorm(1 - T1s)
setbased_corEst1 <- cor(T1s[1:cor.nperm, ])


rm(T1s)


output <- list()
n.perm <- 10000
tmp <- svd(trait.cor)


eigenvalue <- tmp$d
eigenvalue[eigenvalue < max(eigenvalue) / input.cutoff] <- 0
n.eigen = sum(eigenvalue!=0)
# different scenarios

# mu = as.vector(eV$vectors[,1:10] %*% as.matrix(rnorm(10,0,1)))

if(idpgroup == "Freesurfer") {
    delta <- c(0,0.05,0.1, 0.15,0.2, 0.25,0.3,0.35, 0.4,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,0.85, 0.9)
} else {
    delta <- c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65,0.7, 0.8, 0.9, 1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4)
}

for (j in 1:length(delta)) {
    
    T2s <- matrix(NA, n.perm, 11)
    mu = as.vector(tmp$u[,1:n.eigen] %*% as.matrix( eigenvalue[1:n.eigen] * delta[j]))

  
    for (i in 1:n.perm) {
        tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
        tmp.z <- as.vector(tmp.z)
        tmp.z <- tmp.z + mu

        pSum <- Sum_fast(tmp.z, sum_denom)
        pSSU <- SumSqU_fast(tmp.z, beta1, alpha1, d1)
        pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
        pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
        pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
        pAT50 <- MultiXcan_fast(tmp.z, t.mat50, multi.df50)$multi.p

        univ <- 2 * pnorm(-max(abs(tmp.z))) * dim(trait.cor)[1]
        pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p

        x1 <- matrix(tmp.z, ncol = length(tmp.z), nrow = 1)
        T <- W %*% Sigma %*% t(x1)
        T <- (T * T) / (W %*% Sigma %*% t(W))
        SHom <- pchisq(T[1, 1], df = 1, ncp = 0, lower.tail = F)
        if (idpgroup == "Freesurfer_volume") {
            minP <- matz(tmp.z, trait.cor)
        } else {
            minP <- 1 # PowerUniv(tmp.z, trait.cor)
        }

        omni_stat <- min(pAT1, pAT10, pAT30, pAT50)

        omni_p <- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(qnorm(1 - omni_stat), 4), sigma = setbased_corEst1)[1]

        # minP <- PowerUniv(tmp.z, trait.cor)
        T2s[i, ] <- c(pSum, pSSU, univ, omni_p, pAT1, pAT10, pAT30, pAT50, pATall, SHom, minP)
    }

    output[[j]] <- T2s
}

save(output,delta, file = paste("/gpfs/research/chongwu/Chong/Multitrait-snp/simulation/res/PC2_",input.cutoff, "_", idpgroup, ".RData", sep = ""))

# sswrite.table(res, paste("/home/panwei/wuxx0845/Multitrait-snp/simulation/null/", setting,"_",idpgroup, "_", job, ".txt", sep = ""), row.names = F, quote = F)

cat(time)

cutoff <- 0.05

res <- matrix(NA, 8, length(output))
for (i in 1:length(output)) {
    tmp <- output[[i]]
    res[, i] <- colSums(tmp < cutoff) / dim(tmp)[1]
}

rownames(res) <- c("Sum", "SSU", "univ", "aMAT", "MAT(1)", "MAT(10)", "MAT(30)", "MAT(50)")
colnames(res) <- delta
res
