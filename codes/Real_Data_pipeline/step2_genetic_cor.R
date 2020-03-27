setwd("/home/panwei/wuxx0845/Multitrait-snp/ldsc/")

library(data.table)
library(readr)
# library(MTAR)
suppressMessages(library("plink2R"))

suppressMessages(library("optparse"))

args <- commandArgs(TRUE)
job <- (eval(parse(text = args[[1]])))
idpgroup <- (as.character(args[[2]]))




# job <- 5
#idpgroup = "rfMRI"

id <- readRDS("/home/panwei/wuxx0845/UKBiobankImage/IDP_id.rds")

if (idpgroup == "dMRI") {
  id <- id[grepl("IDP_dMRI", id[, 2]),]
}

if (idpgroup == "Freesurfer") {
  id <- id[grepl("volume_", id[, 2]) | grepl("_thickness", id[, 2]) | grepl("_area", id[, 2]),]
  id <- id[8:dim(id)[1],]
}

if (idpgroup == "rfMRI") {
  id <- id[grepl("NODEamps100_", id[, 2]) | grepl("NODEamps25_", id[, 2]),]
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

job.list <- job.list[!is.na(job.list[, 1]),]

# parallel the jobs
each <- ceiling(dim(job.list)[1] / 505)
last.job.num <- floor(dim(job.list)[1] / each)

if (job == last.job.num) {
  job.list <- job.list[(job * each + 1):dim(job.list)[1],]
} else {
  job.list <- job.list[(job * each + 1):(job * each + each),]
}

final.res <- as.data.frame(matrix(NA, dim(job.list)[1], 3))

final.res[, 1:2] <- job.list
for (i in 1:dim(job.list)[1]) {
  #
  tryCatch({
    idp1 <- job.list[i, 1]
    idp2 <- job.list[i, 2]

    idp1name <- gsub(".sumstats", "", idp1)
    idp2name <- gsub(".sumstats", "", idp2)

    output <- paste("./tmpresults/", idpgroup, "_", job, "_", idp1name, "_", idp2name, sep = "")


    idpdir <- "/home/panwei/wuxx0845/Multitrait-snp/ldsc_IDP/"
    command <- paste("./ldsc.py --rg ", idpdir, idp1, ",", idpdir, idp2, "  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ", output, sep = "")

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

file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, "_", job, ".rds", sep = "")
saveRDS(final.res, file.out)


file.out <- paste("/home/panwei/wuxx0845/Multitrait-snp/ldsc/results/res_", idpgroup, "_", job, ".txt", sep = "")
write.table(final.res, file.out, row.names = F, col.names = F, quote = F)

output <- paste("./tmpresults/", idpgroup, "_", job, "_*", sep = "")
system(paste("rm ", output, sep = ""))


idpgroup <- "dMRI"

if (idpgroup == "T1_brain_vol") {
  id <- id[1:10,]
}

if (idpgroup == "T2_star") {
  id <- id[grepl("IDP_SWI_T2star", id[, 2]),]
}

if (idpgroup == "ICA") {
  dim(id[grepl("ICA_", id[, 2]),])
  id <- id[grepl("ICA_", id[, 2]),]
}

if (idpgroup == "dMRI_FA") {
  id <- id[grepl("_FA_", id[, 2]),]
}
if (idpgroup == "dMRI_MD") {
  id <- id[grepl("_MD_", id[, 2]),]
}

if (idpgroup == "dMRI_MO") {
  id <- id[grepl("_MO_", id[, 2]),]
}

if (idpgroup == "dMRI_L1") {
  id <- id[grepl("_L1_", id[, 2]),]
}

if (idpgroup == "dMRI_L2") {
  id <- id[grepl("_L2_", id[, 2]),]
}

if (idpgroup == "dMRI_L3") {
  id <- id[grepl("_L3_", id[, 2]),]
}

if (idpgroup == "dMRI_ICVF") {
  id <- id[grepl("_ICVF_", id[, 2]),]
}

if (idpgroup == "dMRI_ICVF") {
  id <- id[grepl("_ICVF_", id[, 2]),]
}
