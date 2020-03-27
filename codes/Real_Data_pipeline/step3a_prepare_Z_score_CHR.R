setwd("/home/panwei/wuxx0845/UKBiobankImage")


library(data.table)
library(readr)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))

pos = readRDS("processed_hap_pos.rds")

pos$usedindx = pos[,1]==job

pos2 = pos[pos$usedindx,]

id = readRDS("IDP_id.rds")
id = id[grepl("T1_FAST_ROIs",id[,2]),]

id$saveout = gsub(".txt.gz","",id$file)

id$Zpath = paste("Z_",gsub(".txt.gz","",id$file),".rds",sep="")

#test; using the first two 

Zscore = matrix(NA,dim(pos2)[1],dim(id)[1])

for (i in 1:dim(id)[1]) { #dim(id)[1]
    file.name = id[i,"Zpath"]
    file.out = id[i,"saveout"]

    z.tmp = readRDS(file.name)
    z.tmp = z.tmp[pos$usedindx,]
    Zscore[,i] = z.tmp
}

saveRDS(Zscore,paste("/home/panwei/wuxx0845/Multitrait-snp/data/IDP_T1_FAST_ROIs_CHR",job,".rds",sep=""))

