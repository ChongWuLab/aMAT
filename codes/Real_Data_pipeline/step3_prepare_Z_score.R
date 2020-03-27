setwd("/home/panwei/wuxx0845/UKBiobankImage")


library(data.table)
library(readr)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))

#idpgroup = "smallIDPcombined"

pos <- readRDS("processed_hap_pos.rds")

each <- ceiling(dim(pos)[1] / 263)
last.job.num <- floor(dim(pos)[1] / each)

id <- readRDS("IDP_id.rds")

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


if (idpgroup == "smallIDPcombined") {
    id2 <- id[1:25, ]
    id3 <- id[grepl("IDP_T1_FIRST_left", id[, 2]), ]
    id3 <- id3[8:14, ]
    id4 <- id[grepl("IDP_SWI_T2star", id[, 2]), ]
    id5 <- id[grepl("ICA_", id[, 2]), ]
    id6 <- id[grepl("tfMRI_", id[, 2]), ]

    id <- rbind(id2, id3)
    id <- rbind(id,id4)
    id <- rbind(id,id5)
    id <- rbind(id,id6)
}

id$saveout <- gsub(".txt.gz", "", id$file)

id$Zpath <- paste("Z_", gsub(".txt.gz", "", id$file), ".rds", sep = "")

# test; using the first two
for (i in 1:dim(id)[1]) { # dim(id)[1]
    file.name <- id[i, "Zpath"]
    file.out <- id[i, "saveout"]

    z.tmp <- readRDS(file.name)
    for (job in 0:last.job.num) {
        if (job == last.job.num) {
            usedindx <- (job * each + 1):dim(pos)[1]
        } else {
            usedindx <- (job * each + 1):(job * each + each)
        }
        z.tmp2 <- z.tmp[usedindx, ]
        saveRDS(z.tmp2, paste("/home/panwei/wuxx0845/Multitrait-snp/tmp/tmp_ROIs_",idpgroup,"_",i,"_", job, ".rds", sep = ""))
    }
}

for (job in 0:last.job.num) {
    pos$usedindx <- FALSE
    if (job == last.job.num) {
        pos$usedindx[(job * each + 1):dim(pos)[1]] <- TRUE
    } else {
        pos$usedindx[(job * each + 1):(job * each + each)] <- TRUE
    }

    # pos$usedindx = pos[,1]==job

    pos2 <- pos[pos$usedindx, ]

    Zscore <- matrix(NA, dim(pos2)[1], dim(id)[1])
    colnames(Zscore) <- paste("IDP_", id$saveout, sep = "")
    for (i in 1:dim(id)[1]) { # dim(id)[1]
      
        z.tmp <- readRDS( paste("/home/panwei/wuxx0845/Multitrait-snp/tmp/tmp_ROIs_",idpgroup,"_",i,"_", job, ".rds", sep = "") )
        Zscore[, i] <- z.tmp
    }
    Zscore <- cbind(pos2, Zscore)
    saveRDS(Zscore, paste("/home/panwei/wuxx0845/Multitrait-snp/data/",idpgroup,"_split", job, ".rds", sep = ""))
}


output <- paste("home/panwei/wuxx0845/Multitrait-snp/tmp/tmp_ROIs_", idpgroup, "_*", sep = "")
system(paste("rm ", output, sep = ""))