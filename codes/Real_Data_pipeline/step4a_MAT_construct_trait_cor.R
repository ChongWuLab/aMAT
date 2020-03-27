setwd("/home/panwei/wuxx0845/UKBiobankImage")


library(data.table)
library(mvtnorm)
# library(MTAR)
source("/home/panwei/wuxx0845/Multitrait-snp/multitrait_test_support.R")
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))


job = 1
id <- readRDS("IDP_id.rds")

#idpgroup <- "dMRI"

# need to use 504 cores; for others, using 264 cores.
if (idpgroup == "T1_FAST_ROIs") {
    id <- id[grepl("T1_FAST_ROIs", id[, 2]), ]

    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    id <- id[!id$usedIDP=="IDP_0061", ]

    Zscore <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/IDP_T1_FAST_ROIs_split", job, ".rds", sep = ""))
    snp.inf <- Zscore[, 1:10]
    Zscore <- Zscore[, 11:dim(Zscore)[2]]
     Zscore <- Zscore[, id$usedIDP]
    trait.cor <- readRDS("/home/panwei/wuxx0845/Multitrait-snp/data/T1_FAST_ROIs_cor.rds")

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}

dMRI.set <- c("dMRI", "dMRI_ProbtrackX", "dMRI_FA", "dMRI_MD", "dMRI_MO", "dMRI_L1", "dMRI_L2", "dMRI_L3", "dMRI_ICVF", "dMRI_OD", "dMRI_ISOVF")

if (idpgroup %in% dMRI.set) {
    if (idpgroup == "dMRI") {
        id <- id[grepl("IDP_dMRI", id[, 2]), ]
    }

    if (idpgroup == "dMRI_ProbtrackX") {
        id <- id[grepl("dMRI_ProbtrackX", id[, 2]), ]
    }

    if (idpgroup == "dMRI_FA") {
        id <- id[grepl("_FA_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_MD") {
        id <- id[grepl("_MD_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_MO") {
        id <- id[grepl("_MO_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_L1") {
        id <- id[grepl("_L1_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_L2") {
        id <- id[grepl("_L2_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_L3") {
        id <- id[grepl("_L3_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_ICVF") {
        id <- id[grepl("_ICVF_", id[, 2]), ]
    }

    if (idpgroup == "dMRI_OD") {
        id <- id[grepl("_OD_", id[, 2]), ]
    }
    if (idpgroup == "dMRI_ISOVF") {
        id <- id[grepl("_ISOVF_", id[, 2]), ]
    }
    # start here tomorrow
    remain <- readRDS("/home/panwei/wuxx0845/Multitrait-snp/data/remain_dMRI_cor.rds")
    remain.tmp <- c(remain[, 1], remain[, 2])
    remain.tmp <- table(remain.tmp)
    remain.tmp <- sort(remain.tmp)
    remove.idp <- names(remain.tmp)[remain.tmp > 600]
    remove.idp <- gsub(".sumstats", "", remove.idp)

    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    id <- id[!id$usedIDP %in% remove.idp, ]


    Zscore <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/dMRI_split", job, ".rds", sep = ""))

    snp.inf <- Zscore[, 1:10]
    Zscore <- Zscore[, 11:dim(Zscore)[2]]
    Zscore <- Zscore[, id$usedIDP]

    trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "dMRI", "_cor.rds", sep = ""))

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}

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

    remain <- readRDS("/home/panwei/wuxx0845/Multitrait-snp/data/remain_Freesurfer_cor.rds")
    remain.tmp <- c(remain[, 1], remain[, 2])
    remain.tmp <- table(remain.tmp)
    remain.tmp <- sort(remain.tmp)
    remove.idp <- names(remain.tmp)[remain.tmp > 470]
    remove.idp <- gsub(".sumstats", "", remove.idp)

    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    id <- id[!id$usedIDP %in% remove.idp, ]


    Zscore <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/Freesurfer_split", job, ".rds", sep = ""))

    snp.inf <- Zscore[, 1:10]
    Zscore <- Zscore[, 11:dim(Zscore)[2]]
    Zscore <- Zscore[, id$usedIDP]

    trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "Freesurfer", "_cor.rds", sep = ""))

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}

rfMRI.set <- c("rfMRI", "rfMRI_25", "rfMRI_100")
if (idpgroup %in% rfMRI.set) {
    if (idpgroup == "rfMRI") {
        id <- id[grepl("NODEamps100_", id[, 2]) | grepl("NODEamps25_", id[, 2]), ]
    }

    if (idpgroup == "rfMRI_25") {
        id <- id[ grepl("NODEamps25_", id[, 2]), ]
    }

    if (idpgroup == "rfMRI_100") {
        id <- id[grepl("NODEamps100_", id[, 2]), ]
    }



    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    Zscore <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/rfMRI_split", job, ".rds", sep = ""))

    snp.inf <- Zscore[, 1:10]
    Zscore <- Zscore[, 11:dim(Zscore)[2]]
    Zscore <- Zscore[, id$usedIDP]

    trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "rfMRI", "_cor.rds", sep = ""))

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}


smallIDP.set <- c("T1_brain_subcortical", "T2_star", "ICA", "tfMRI", "T1_brain_vol", "T1_Subcortical", "T1_Subcortical_L_plus_R", "T2_star_L_plus_R", "T2_star_combined")
if (idpgroup %in% smallIDP.set) {
    if (idpgroup == "T1_brain_vol") {
        id <- id[1:10, ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T1_brain_subcortical", "_cor.rds", sep = ""))
    }

    if (idpgroup == "T1_Subcortical") {
        id <- id[11:25, ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T1_brain_subcortical", "_cor.rds", sep = ""))
    }

    if (idpgroup == "T1_Subcortical_L_plus_R") {
        id <- id[grepl("IDP_T1_FIRST_left", id[, 2]), ]
        id <- id[8:14, ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T1_brain_subcortical", "_cor.rds", sep = ""))
    }

    if (idpgroup == "T1_brain_subcortical") {
        id2 <- id[1:25, ]
        id3 <- id[grepl("IDP_T1_FIRST_left", id[, 2]), ]
        id3 <- id3[8:14, ]
        id <- rbind(id2, id3)
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T1_brain_subcortical", "_cor.rds", sep = ""))
    }

    if (idpgroup == "T2_star_combined") {
        id <- id[grepl("IDP_SWI_T2star", id[, 2]), ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T2_star", "_cor.rds", sep = ""))
    }


    if (idpgroup == "T2_star_L_plus_R") {
        id <- id[grepl("IDP_SWI_T2star", id[, 2]), ]
        id <- id[15:21, ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T2_star", "_cor.rds", sep = ""))
    }

    if (idpgroup == "T2_star") {
        id <- id[grepl("IDP_SWI_T2star", id[, 2]), ]
        id <- id[1:14, ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "T2_star", "_cor.rds", sep = ""))
    }

    if (idpgroup == "ICA") {
        id <- id[grepl("ICA_", id[, 2]), ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "ICA", "_cor.rds", sep = ""))
    }

    if (idpgroup == "tfMRI") {
        id <- id[grepl("tfMRI_", id[, 2]), ]
        trait.cor <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/", "tfMRI", "_cor.rds", sep = ""))
    }

    remain <- readRDS("/home/panwei/wuxx0845/Multitrait-snp/data/remain_tfMRI_cor.rds")
    remain.tmp <- c(remain[, 1], remain[, 2])
    remain.tmp <- table(remain.tmp)
    remain.tmp <- sort(remain.tmp)
    remove.idp <- names(remain.tmp)[remain.tmp > 10]
    remove.idp <- gsub(".sumstats", "", remove.idp)

    id$saveout <- gsub(".txt.gz", "", id$file)
    id$usedIDP <- paste("IDP_", id$saveout, sep = "")
    id <- id[!id$usedIDP %in% remove.idp, ]

    Zscore <- readRDS(paste("/home/panwei/wuxx0845/Multitrait-snp/data/smallIDPcombined_split", job, ".rds", sep = ""))

    snp.inf <- Zscore[, 1:10]
    Zscore <- Zscore[, 11:dim(Zscore)[2]]
    Zscore <- Zscore[, id$usedIDP]

    colnames(trait.cor) <- gsub(".sumstats", "", colnames(trait.cor))
    rownames(trait.cor) <- gsub(".sumstats", "", rownames(trait.cor))

    trait.cor <- trait.cor[id$usedIDP, ]
    trait.cor <- trait.cor[, id$usedIDP]
}


trait.diag <- readRDS("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_genetic_cor_diagonal.rds")
trait.diag[, 3] <- gsub(".sumstats", "", trait.diag[, 3])
rownames(trait.diag) <- trait.diag[, 3]

trait.diag <- trait.diag[rownames(trait.cor), ]
trait.diag[is.na(trait.diag[, 4]), 4] <- 1
diag(trait.cor) <- trait.diag[, 4]


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
eigenvalue[eigenvalue < 0] <- 0
eigenvalue[eigenvalue != 0] <- 1 / eigenvalue[eigenvalue != 0]
multi.dfall <- sum(eigenvalue != 0)
D <- diag(eigenvalue)
t.matall <- tmp$u %*% D %*% t(tmp$v)
invhalf.matall <- tmp$u %*% sqrt(D) %*% t(tmp$v)


#methods.cor <- omnibus_cor2(  trait.cor, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall,n.perm = 10000)

methods.cor2 <- omnibus_cor3(
    trait.cor, sum_denom, beta1, alpha1, d1, t.mat1, multi.df1, t.mat10, multi.df10, t.mat30, multi.df30, t.mat50, multi.df50, t.matall, multi.dfall,
    n.perm = 10000
)

saveRDS(methods.cor2, paste("/home/panwei/wuxx0845/Multitrait-snp/MAT_cor/cor_", idpgroup, ".rds", sep = ""))
