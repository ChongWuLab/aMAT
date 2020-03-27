setwd("/home/panwei/wuxx0845/UKBiobankImage/")


library(data.table)
library(R.utils)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))

id.list <- list.files("/home/panwei/public/Oxford_BIG/GWAS_datasets_3144/", recursive = FALSE)


each <- ceiling(length(id.list) / 70)
last.job.num <- floor(length(id.list) / each)

if (job == last.job.num) {
    id.list <- id.list[(job * each + 1):length(id.list)]
} else {
    id.list <- id.list[(job * each + 1):(job * each + each)]
}

pos <- fread("/home/panwei/wuxx0845/UKBiobankImage/positions.txt")
pos <- as.data.frame(pos)
keep <- nchar(pos$REF) == 1 & nchar(pos$ALT) == 1

a1 <- pos$REF
a2 <- pos$ALT
keep2 <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
keep <- keep & keep2

for (indx in id.list) {
    #indx <- id.list[2]
    #indx = "0045.txt.gz"
    output.txt <- gsub(".gz", "", indx)
    output <- gsub(".txt", "", output.txt)
    file.in <- paste("/home/panwei/public/Oxford_BIG/GWAS_datasets_3144/", indx, sep = "")
    file.out <- paste("/home/panwei/wuxx0845/UKBiobankImage/", output.txt, sep = "")

    gunzip(file.in, file.out, temporary = F, remove = F, overwrite = T)

    tmp <- fread(file.out)
    tmp <- as.data.frame(tmp)
    tmp$Z <- tmp$BETA/tmp$SEBETA
    tmp = tmp[keep,]
    tmp2 = as.matrix(tmp$Z)
    saveRDS(tmp2,paste("Z_",output,".rds",sep=""))
    system(paste("rm ",file.out, sep=""))
}
