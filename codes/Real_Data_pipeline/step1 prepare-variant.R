setwd("/home/chong/data/UKBiobank")

library(data.table)
suppressMessages(library("plink2R"))
suppressMessages(library("optparse"))


pos <- fread("positions.txt")
pos <- as.data.frame(pos)
keep <- nchar(pos$REF) == 1 & nchar(pos$ALT) == 1

a1 <- pos$REF
a2 <- pos$ALT
keep2 <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
keep <- keep & keep2

pos = pos[keep,]
#pos$id =paste(pos$CHROM,":",pos$POS,sep="")
pos$id2 =paste(pos$CHROM,":",pos$POS,":",pos$REF,":",pos$ALT,sep="")

variant = fread("variants.tsv")
variant= as.data.frame(variant)

variant = variant[,c(1:6,10)]
#variant$id = paste(variant$chr,":",variant$pos,sep="")

#sum(pos$id %in% variant$id)
sum(pos$id2 %in% variant[,1])

pos$rsid = NA

rownames(pos) = pos$id2
rownames(variant) = variant[,1]

pos$indx = pos$id2 %in% variant[,1]

tmp.name = pos[pos$indx,"id2"]

variant2= variant[tmp.name,]

pos[pos$indx,"rsid"] = variant2[,"rsid"]

saveRDS(pos,"tmp_pos.rds")


pos = readRDS("tmp_pos.rds")
length(unique(pos$id))

hapmapsnp = fread("w_hm3.snplist")
hapmapsnp = as.data.frame(hapmapsnp)

sum(pos$rsid %in% hapmapsnp[,1])

pos$hapmapindx = pos$rsid %in% hapmapsnp[,1]
pos$haptmp = !duplicated(pos$rsid)
pos$hapmapindx = pos$haptmp & pos$hapmapindx

sum(pos$hapmapindx)
saveRDS(pos,"processed_hap_pos.rds")

big.sum = fread("BIG_summary_stats_files.csv")
big.sum = as.data.frame(big.sum)

big.sum = big.sum[,c(1,2,4)]

hert = fread("hertibability.csv")
hert = as.data.frame(hert)

rownames(big.sum) = big.sum[,2]
rownames(hert) = hert[,1]


hert =hert[rownames(big.sum),]

big = cbind(big.sum,hert)


saveRDS(big,"IDP_id.rds")

big2 = big[2655:2715,]
big2 = big2[2:(dim(big2)[1]-2),]
sum(big.sum[,2] %in% "Freesurfer_volume")
grepl("T1_FAST_ROIs",big[,2])








tar -jxvf eur_w_ld_chr.tar.bz2
unzip -o pgc.cross.bip.zip
unzip -o pgc.cross.scz.zip
bunzip2 w_hm3.snplist.bz2
munge_sumstats.py --sumstats pgc.cross.SCZ17.2013-05.txt --N 17115 --out scz --merge-alleles w_hm3.snplist

./munge_sumstats.py --sumstats pgc.cross.BIP11.2013-05.txt --N 11810 --out bip --merge-alleles w_hm3.snplist

# LD Score Regression
./ldsc.py --rg scz.sumstats2.gz,bip.sumstats2.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_bip

./ldsc.py --rg scz.sumstats,bip.sumstats --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_bip


less scz_bip.log
