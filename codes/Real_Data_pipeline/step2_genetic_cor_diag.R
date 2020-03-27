setwd("/home/panwei/wuxx0845/Multitrait-snp/ldsc/")

library(data.table)
library(readr)
# library(MTAR)
suppressMessages(library("plink2R"))

suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))


id <- readRDS("/home/panwei/wuxx0845/UKBiobankImage/IDP_id.rds")
id$ID <- gsub(".txt.gz","",id$file)

idp = list.files(path = "/home/panwei/wuxx0845/Multitrait-snp/ldsc_IDP/")
idp = idp[1:1452]

used.id = gsub(".sumstats","",idp)
used.id = gsub("IDP_","",used.id)

id = id[id$ID %in% used.id,]
id = id[,c("IDP","ID")]
rownames(id) = id$ID
id = id[used.id,]

id$input = idp
final.res = as.data.frame(matrix(NA,dim(id)[1],4))
final.res[,1:3] = id
for (i in 1:dim(id)[1]) { #
    tryCatch({
        idp1 <- id[i, 3]

        idp1name <- id[i,2]

        output <- paste("./tmpresults/", idp1name,  sep = "")
        idpdir <- "/home/panwei/wuxx0845/Multitrait-snp/ldsc_IDP/"
        command <- paste("./ldsc.py --h2 ", idpdir, idp1, "  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ", output, sep = "")

        system(command)

        tmp <- read_file(paste(output, ".log", sep = ""))

        tmp <- gsub("\n\nRatio.*", "", tmp)
        tmp <- gsub(".*\nIntercept: ", "", tmp)
        tmp <- gsub(" \\(.*", "", tmp)
        tmp <- as.numeric(tmp)

        final.res[i, 4] <- tmp
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

file.out <- "/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_genetic_cor_diagonal.rds"
saveRDS(final.res, file.out)
