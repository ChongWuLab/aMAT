
library(xtable)

library(data.table)
read.file <- function(file.name) {
    file <- try(fread(file.name))
    if (class(file) == "try-error") {
        # cat("Caught an error during fread, trying read.table.\n")
        cat(file.name)
    }
    file
}

read.file2 <- function(file.name) {
    file <- try(readRDS(file.name))
    if (class(file) == "try-error") {
        # cat("Caught an error during fread, trying read.table.\n")
        cat(file.name)
    }
    file
}

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))

idpgroup <- "T1_FAST_ROIs"
job = 262
setwd("/home/panwei/wuxx0845/Multitrait-snp/step4Res")
output <- NULL
# time = matrix(0,,5)

out.i <- 0
seed <- NULL
for (index in 0:job)
{
    #file.name <- paste("res_", idpgroup,"_",index, ".txt", sep = "")

    #res.temp <- read.file(file.name)
    file.name <- paste("sig_res_", idpgroup,"_",index, ".rds", sep = "")
    res.temp <- read.file2(file.name)
    if (is.data.frame(res.temp)) {
        output <- rbind(output,res.temp)
        out.i <- out.i + 1
        if(out.i%%10 == 0) {
            cat(out.i, " ")

        }
    } else {
        seed <- c(seed, index)
        cat("missing: ",index,"\n")
    }
}

write.table(output, paste("/home/panwei/wuxx0845/Multitrait-snp/resAna/res2_", idpgroup, ".txt", sep = ""), row.names = F, quote = F)
