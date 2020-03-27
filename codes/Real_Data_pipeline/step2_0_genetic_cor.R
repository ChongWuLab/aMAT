setwd("/home/panwei/wuxx0845/UKBiobankImage")

library(data.table)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))



ldscore <- NULL
for (i in 1:22) {
    tmp <- fread(paste("/home/panwei/wuxx0845/UKBiobankImage/eur_w_ld_chr/", i, ".l2.ldscore", sep = ""))
    tmp <- as.data.frame(tmp)
    ldscore <- rbind(ldscore, tmp)
}


pos <- readRDS("processed_hap_pos.rds")
pos$ldsnp <- pos$rsid %in% ldscore[, "SNP"]
pos$usedindx <- pos$hapmapindx & pos$ldsnp

pos2 <- pos[pos$usedindx, ]
pos2 <- pos2[, c(6, 3, 4)]

id <- readRDS("IDP_id.rds")
#indx <- c(1:869, 2642:3143)

indx <- c(870:950)
id <- id[indx, ]

id$saveout <- gsub(".txt.gz", "", id$file)

id$Zpath <- paste("Z_", gsub(".txt.gz", "", id$file), ".rds", sep = "")

ldscore.used <- ldscore[ldscore[, "SNP"] %in% pos2[, 1], ]


for (i in 1:dim(id)[1]) {
    tryCatch({
        file.name <- id[i, "Zpath"]
        file.out <- id[i, "saveout"]

        z.tmp <- readRDS(file.name)
        z.tmp <- z.tmp[pos$usedindx, ]
        # z.tmp =  formatC(z.tmp, digits = 8, format = "f")
        tmp <- cbind(pos2, z.tmp)
        tmp$N <- 8411
        colnames(tmp) <- c("SNP", "A1", "A2", "Z", "N")

        write.table(tmp, paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc_IDP/IDP_", file.out, ".sumstats", sep = ""), row.names = F, quote = F)
    }, error = function(e) {
        cat("job id: ", i, "   ")
        cat("ERROR :", conditionMessage(e), "\n")
    })
}