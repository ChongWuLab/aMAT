#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

setwd("/gpfs/research/chongwu/Chong/UKBiobankImage")
library("MASS")
library("Matrix")
require(compiler)
enableJIT(4)

library(data.table)
library(mvtnorm)
# library(MTAR)
source("/gpfs/research/chongwu/Chong/Multitrait-snp/multitrait_test_support.R")
source("/gpfs/research/chongwu/Chong/Multitrait-snp/FunctionSet.R")

suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))


# pos <- readRDS("processed_hap_pos.rds")

# job

id <- readRDS("IDP_id.rds")
idpgroup <- "Freesurfer"

Freesurfer.set <- c("Freesurfer", "Freesurfer_volume", "Freesurfer_area", "Freesurfer_thickness")
if (idpgroup %in% Freesurfer.set) {
    if (idpgroup == "Freesurfer") {
        id <- id[grepl("volume_", id[, 2]) | grepl("_thickness", id[, 2]) | grepl("_area", id[, 2]), ]
        id <- id[8:dim(id)[1], ]
    }

    if (idpgroup == "Freesurfer_volume") {
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

#trait.cor2 = trait.cor + diag(0.01,dim(trait.cor)[1])
#tmp <- svd(trait.cor2)
#tmp$d

tmp <- svd(trait.cor)
eigenvalue <- tmp$d
#if (setting == "setting1") {
#    eigenvalue[eigenvalue < 0.005] <- 1e-6
#}

#eigenvalue[eigenvalue < 0.01] <- 1e-6

#if (setting == "setting2") {
#    eigenvalue[eigenvalue < 0.01] <- 1e-6
#}

#D <- diag(eigenvalue)
#trait.cor <- tmp$u %*% D %*% t(tmp$v)

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


# methods.cor <- omnibus_cor2(  trait.cor, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall,n.perm = 10000)

SampleSize <- rep(8411, dim(trait.cor)[1])

Wi <- matrix(SampleSize, nrow = 1)
sumW <- sqrt(sum(Wi^2))
W <- Wi / sumW
Sigma <- ginv(trait.cor)
n.perm <- 10000


tmp <- svd(trait.cor)
eigenvalue <- tmp$d
D <- diag(eigenvalue)
CovSsqrt <- tmp$u %*% sqrt(D) %*% t(tmp$v)

#eV <- eigen(trait.cor)
#eV$values[eV$values < 0] <- 0


T1s <- matrix(NA, n.perm, 4)
start.time <- proc.time()[3]

for (i in 1:n.perm) {
    tmp.z <- rnorm(dim(CovSsqrt)[1], 0, 1) %*% CovSsqrt
    tmp.z <- as.vector(tmp.z)
    pAT1 <- MultiXcan_fast(tmp.z, t.mat1, multi.df1)$multi.p
    pAT10 <- MultiXcan_fast(tmp.z, t.mat10, multi.df10)$multi.p
    pAT30 <- MultiXcan_fast(tmp.z, t.mat30, multi.df30)$multi.p
    pAT50 <- MultiXcan_fast(tmp.z, t.mat50, multi.df50)$multi.p
    #pATrev <- MultiXcan_fast(tmp.z, t.matallrev, multi.dfallrev)$multi.p

    T1s[i, ] <- c(pAT1, pAT10, pAT30, pAT50)
}

colSums(T1s < 0.01,na.rm=T) / n.perm
T1s.save <- T1s
T1s[T1s > 0.99] <- 0.99
T1s <- qnorm(1 - T1s)
setbased_corEst1 <- cor(T1s)


n.perm = 1000000
res <- matrix(NA, n.perm, 9)

start.time = proc.time()[3]
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
    #pATrev <- MultiXcan_fast(tmp.z, t.matallrev, multi.dfallrev)$multi.p
    pATall <- MultiXcan_fast(tmp.z, t.matall, multi.dfall)$multi.p

    omni_stat <- min(pAT1, pAT10, pAT30, pAT50)
    aMAT<- 1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(qnorm(1 - omni_stat), 4), sigma = setbased_corEst1)[1]

    x1 <- matrix(tmp.z, ncol = length(tmp.z), nrow = 1)
    T <- W %*% Sigma %*% t(x1)
    T <- (T * T) / (W %*% Sigma %*% t(W))
    SHom <- pchisq(T[1, 1], df = 1, ncp = 0, lower.tail = F)

    res[i, ] <- c(pSum, pSSU,pATall, SHom, pAT1, pAT10, pAT30, pAT50, aMAT)
}
proc.time()[3] - start.time


final.out = matrix(NA,11,9)

cut.off = c(0.05,0.01,5e-3,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8)
rownames(final.out) = cut.off
colnames(final.out) = c("Sum","SSU","chi-square", "Hom", "MAT1", "MAT10", "MAT30", "MAT50", "aMAT")

for(i in 1:length(cut.off)) {
    final.out[i,] = colSums(res < cut.off[i],na.rm=T)
}

n.sim = colSums(!is.na(res))

save(final.out, n.sim,file=paste("/gpfs/research/chongwu/Chong/Multitrait-snp/simulation/res/null_Freesurf/null_", job,"_",job.id, ".RData", sep = ""))

