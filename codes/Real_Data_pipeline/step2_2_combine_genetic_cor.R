setwd("/home/panwei/wuxx0845/Multitrait-snp/ldsc/")

library(data.table)
library(readr)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))



# job <- 5

get.cor <- function(idpgroup) {
   # idpgroup <- "rfMRI"
    id <- readRDS("/home/panwei/wuxx0845/UKBiobankImage/IDP_id.rds")

    # idpgroup <- "T1_brain_subcortical"
    if (idpgroup == "dMRI") {
        id <- id[grepl("IDP_dMRI", id[, 2]), ]
    }

    if (idpgroup == "Freesurfer") {
        id <- id[grepl("volume_", id[, 2]) | grepl("_thickness", id[, 2]) | grepl("_area", id[, 2]), ]
        id <- id[8:dim(id)[1], ]
    }

    if (idpgroup == "rfMRI") {
        id <- id[grepl("NODEamps100_", id[, 2]) | grepl("NODEamps25_", id[, 2]), ]
    }


    if (idpgroup == "T1_brain_subcortical") {
        id2 <- id[1:25, ]
        id3 <- id[grepl("IDP_T1_FIRST_left", id[, 2]), ]
        id3 <- id3[8:14, ]
        id <- rbind(id2, id3)
    }

    if (idpgroup == "T2_star") {
        id <- id[grepl("IDP_SWI_T2star", id[, 2]), ]
    }

    if (idpgroup == "ICA") {
        dim(id[grepl("ICA_", id[, 2]), ])
        id <- id[grepl("ICA_", id[, 2]), ]
    }

    if (idpgroup == "tfMRI") {
        dim(id[grepl("tfMRI_", id[, 2]), ])
        id <- id[grepl("tfMRI_", id[, 2]), ]
    }


    id$saveout <- gsub(".txt.gz", "", id$file)

    id$IDPpath <- paste("IDP_", id$saveout, ".sumstats", sep = "")

    id <- id$IDPpath

    cor <- matrix(NA, length(id), length(id))
    colnames(cor) <- rownames(cor) <- id

    id.indx <- cbind(id, 1:length(id))
    rownames(id.indx) <- id
    remain <- NULL

    if (idpgroup != "dMRI" & idpgroup != "Freesurfer" & idpgroup != "rfMRI") {
        file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, ".rds", sep = "")
        tmp <- readRDS(file.out)
        remain <- rbind(remain, tmp)

        indx <- tmp[, 1:2]

        indx$row <- id.indx[tmp[, 1], 2]
        indx$col <- id.indx[tmp[, 2], 2]

        indx <- indx[, 3:4]
        indx <- as.matrix(indx)
        indx <- apply(indx, 2, as.numeric)
        cor[indx] <- tmp[, 3]
        indx2 <- indx[, c(2, 1)]
        cor[indx2] <- tmp[, 3]
        diag(cor) <- 1
    } else {
        if (idpgroup == "dMRI" | idpgroup == "Freesurfer") {
            job.indx <- 0:504
        } else {
            job.indx <- 0:475
        }
        for (job in job.indx) {
            file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, "_", job, ".rds", sep = "")
            tmp <- readRDS(file.out)
            if(job == 475 &idpgroup == "rfMRI") {
                tmp = tmp[!is.na(tmp[,3]),]
            }
            remain <- rbind(remain, tmp)

            indx <- tmp[, 1:2]

            indx$row <- id.indx[tmp[, 1], 2]
            indx$col <- id.indx[tmp[, 2], 2]

            indx <- indx[, 3:4]
            indx <- as.matrix(indx)
            indx <- apply(indx, 2, as.numeric)
            cor[indx] <- tmp[, 3]
            if ( (job == 504 & idpgroup == "Freesurfer") | (job == 475 & idpgroup == "rfMRI") ) {
                indx2 <- indx[c(2, 1)]
            } else {
                indx2 <- indx[, c(2, 1)]
            }
            cor[indx2] <- tmp[, 3]
        }
        diag(cor) <- 1
    }
    
    cat(idpgroup, " number of na ", sum(is.na(cor)), "\n")

    remain = remain[is.na(remain[,3]) & !is.na(remain[,2]),]
    if(dim(remain)[1] > 0 ){
        saveRDS(remain,paste("/home/panwei/wuxx0845/Multitrait-snp/data/remain_", idpgroup, "_cor.rds", sep = ""))
    }

    #saveRDS(cor, paste("/home/panwei/wuxx0845/Multitrait-snp/data/", idpgroup, "_cor.rds", sep = ""))

    # tmp <- svd(cor)
}
get.cor("dMRI")
get.cor("Freesurfer")
get.cor("rfMRI")

get.cor("T1_brain_subcortical")
get.cor("T2_star")
get.cor("ICA")
get.cor("tfMRI")









remain <- remain[is.na(remain[, 3]), ]

for (i in 1:dim(remain)[1]) {
    tmp1 <- remain[i, 1]
    tmp1 <- fread(tmp1)

    tmp1 <- as.data.frame(tmp1)
    tmp2 <- remain[i, 2]
    tmp2 <- fread(tmp2)
    tmp2 <- as.data.frame(tmp2)

    remain[i, 3] <- cor(tmp1[, 4], tmp2[, 4])
}

file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", -1, ".rds", sep = "")
saveRDS(remain, file.out)
