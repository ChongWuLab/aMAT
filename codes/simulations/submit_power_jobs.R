#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)


gamma = c(1,10,20,30,40,50)
setting = c("Freesurfer_area","Freesurfer_thickness","fiveIDPs","twentyfiveIDPs")

j = ceiling(job.id/6)
i = job.id - (j-1) *6

command = paste("R CMD BATCH --no-save --no-restore '--args  job=",gamma[i]," ", setting[j]," ' /gpfs/research/chongwu/Chong/Multitrait-snp/simulation/sim3_power_PC2.R  /gpfs/research/chongwu/Chong/Multitrait-snp/simulation/log/sim3_power_",gamma[i],"_",setting[j],".txt",sep="")

system(command)

