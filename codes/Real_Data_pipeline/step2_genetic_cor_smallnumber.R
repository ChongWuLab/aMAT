setwd("/home/panwei/wuxx0845/Multitrait-snp/ldsc/")

library(data.table)
library(readr)
#library(MTAR)

#suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))

#idpgroup ="ICA"
# job <- 5
id <- readRDS("/home/panwei/wuxx0845/UKBiobankImage/IDP_id.rds")

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

# id <- id[grepl("T1_FAST_ROIs", id[, 2]), ]

id$saveout <- gsub(".txt.gz", "", id$file)

id$IDPpath <- paste("IDP_", id$saveout, ".sumstats", sep = "")

tmp <- id$IDPpath

job.list <- as.data.frame(matrix(NA, length(tmp) * length(tmp) / 2 + 1000, 2))
index <- 0
for (i in 1:(length(tmp) - 1)) {
    index.tmp <- (i + 1):length(tmp)
    index.tmp2 <- index.tmp + index

    job.list[index.tmp2, 1] <- tmp[i]
    job.list[index.tmp2, 2] <- tmp[index.tmp]
    index <- index + length(index.tmp)
}

job.list <- job.list[!is.na(job.list[, 1]), ]

# parallel the jobs


final.res <- as.data.frame(matrix(NA, dim(job.list)[1], 3))

final.res[, 1:2] <- job.list
for (i in 1:dim(job.list)[1]) { #
    tryCatch({
        idp1 <- job.list[i, 1]
        idp2 <- job.list[i, 2]

        idp1name <- gsub(".sumstats", "", idp1)
        idp2name <- gsub(".sumstats", "", idp2)

        output <- paste("./tmpresults/", idpgroup, "_", idp1name, "_", idp2name, sep = "")
        idpdir = "/home/panwei/wuxx0845/Multitrait-snp/ldsc_IDP/"
        command <- paste("./ldsc.py --rg ", idpdir,idp1, ",", idpdir,idp2, "  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ", output, sep = "")

        system(command)

        tmp <- read_file(paste(output, ".log", sep = ""))

        tmp <- gsub("\n\nGenetic Correlation.*", "", tmp)
        tmp <- gsub(".*Total Observed scale gencov", "", tmp)
        tmp <- gsub(".*\nIntercept: ", "", tmp)
        tmp <- gsub(" \\(.*", "", tmp)
        tmp <- as.numeric(tmp)

        final.res[i, 3] <- tmp
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, ".rds", sep = "")
saveRDS(final.res, file.out)


file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, ".txt", sep = "")
write.table(final.res, file.out, row.names = F, col.names = F, quote = F)

output <- paste("./tmpresults/", idpgroup, "_*", sep = "")
system(paste("rm ", output, sep = ""))



