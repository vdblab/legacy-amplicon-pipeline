#!/usr/bin/env Rscript#
#args = commandArgs(trailingOnly=TRUE)

library(dada2); packageVersion("dada2")
library(ShortRead); #To subsample fastq sequences.
#forward
                                        #path <- args[1]
if (exists("snakemake")){
                                        #print(str(snakemake))
    print("parsing snakemake args")
    path <- snakemake@params[["input_dir"]]
    # see https://stackoverflow.com/questions/64101921/
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    fnFs <- unlist(snakemake@input[["readsf"]])
    fnRs <- unlist(snakemake@input[["readsr"]])
    filter_trunclen = c(snakemake@params[["trunclen"]], snakemake@params[["trunclen"]])
    min_asv_len = snakemake@params[["minasvlen"]]
    outpath <- snakemake@params[["outdir_name"]]
    dir.create(outpath)
    dir.create(file.path(outpath, "filtered"))
    ncores <- snakemake@params[["ncores"]]

} else{
    args = commandArgs(trailingOnly=TRUE)
    path <- args[[1]]
    fnFs <- Sys.glob(file.path(path, "*_R1*.fastq.gz"))
    outpath <- "dada2out"
    fnRs <- Sys.glob(file.path(path, "*_R2*.fastq.gz"))
    filter_trunclen = c(180,180) #c(180,180) is default
    min_asv_len = 300
    ncores <- 8
}
#remove unsigned
fnFs <- fnFs[lapply(fnFs,function(x) length(grep("_unsigned",x,value=FALSE))) == 0]
fnRs <- fnRs[lapply(fnRs,function(x) length(grep("_unsigned",x,value=FALSE))) == 0]



asv_len=seq(min_asv_len, 2*filter_trunclen[1]) #c(300,360) is default
print(min(asv_len))
print(max(asv_len))
# filter_trunclen = c(145,151); #c(180,180) is default
# asv_len=seq(240,260); #c(300,360) is default




forwardfileinfo <- file.info(fnFs)
forwardsizes <- forwardfileinfo[,"size"]
reversefileinfo <- file.info(fnRs)
reversesizes <- reversefileinfo[,"size"]#forwardfileinfo[,"size"]


fnFs <- fnFs[which(forwardsizes>10000 | reversesizes>10000)]
fnRs <- fnRs[which(forwardsizes>10000 | reversesizes>10000)]
#get sample names (update depending on filename conventions)
#sample.names <- sapply(strsplit(sapply(strsplit(basename(fnFs), ".fastq.gz"), `[`, 1), "fastq_"), `[`,2)
#sample.names <- sapply(strsplit(sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1), "fastq_"), `[`,2)
#sample.names = sapply(strsplit(basename(fnFs), "_R1.fastq.gz"), `[`, 1);
sample.names = basename(fnFs)
print('sample names')
print(sample.names)

#temporary files for storing filtered data
filtFs <- file.path(outpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(outpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

multithread_var=ncores
randomize_var=F

sprintf("%s - before filterAndTrim Fs",Sys.time())

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=filter_trunclen,
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=multithread_var)

#Capping sampleids up to 50k reads.
sprintf("%s - before capping",Sys.time())

Nread_cap=1e5

file.remove(dir( dirname(filtFs[1]),pattern="temp",full.names=T ))
do_capping=T;#F;

if(do_capping){

  for(i in 1:length(filtFs) ){
    R1_fq_file=filtFs[i];
    R2_fq_file=filtRs[i];

    sampler_R1 <- FastqSampler(R1_fq_file,Nread_cap);
    sampler_R2 <- FastqSampler(R2_fq_file,Nread_cap);

    set.seed(123); sampler_R1_y <- yield(sampler_R1);
    set.seed(123); sampler_R2_y <- yield(sampler_R2);

    temp_R1_fq_file=gsub("fastq.gz","fastq.gz.temp",R1_fq_file);
    temp_R2_fq_file=gsub("fastq.gz","fastq.gz.temp",R2_fq_file);

    writeFastq(sampler_R1_y, temp_R1_fq_file);
    writeFastq(sampler_R2_y, temp_R2_fq_file);

    file.rename(temp_R1_fq_file, R1_fq_file);
    file.rename(temp_R2_fq_file, R2_fq_file);

    close(sampler_R1);
    close(sampler_R2);

  }

}

sprintf("%s - before learnErrors Fs",Sys.time())
errF <- learnErrors(filtFs, multithread=multithread_var, randomize=randomize_var);#imultithread=multithread_var)

sprintf("%s - before learnErrors Rs",Sys.time())
errR <- learnErrors(filtRs, multithread=multithread_var, randomize=randomize_var);#TRUE)

sprintf("%s - before derep F",Sys.time())
derepFs <- derepFastq(filtFs, verbose=TRUE)

sprintf("%s - before derep R",Sys.time())
derepRs <- derepFastq(filtRs, verbose=TRUE)

#names(derepFs) <- sample.names
#names(derepRs) <- sample.names

sprintf("%s - before dada F",Sys.time())
dadaFs <- dada(derepFs, err=errF,multithread=multithread_var);#TRUE)

sprintf("%s - before dada R",Sys.time())
dadaRs <- dada(derepRs, err=errR,multithread=multithread_var);#TRUE)

sprintf("%s - before merge",Sys.time())
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs);#, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

sprintf("%s - before removebimera",Sys.time())

seqtab_prebimera_file = sprintf("%s/seqtab_prebimera.rds",outpath);
saveRDS(seqtab, seqtab_prebimera_file);
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=multithread_var)#, verbose=TRUE)

sprintf("%s - before subsetting by length",Sys.time())
print(str(seqtab.nochim))
print("ASV sequence lengths:")
print(summary(nchar(colnames(seqtab.nochim))))
seqtabfinal <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% asv_len, drop=F]
print(str(seqtabfinal))
#rownames(seqtabfinal) = names(derepFs);
rownames(seqtabfinal) = sample.names; # I changed names to sample.names to avoid error when dealing with 1 dimension data.

sprintf("%s - before writing out noprior files",Sys.time())
outfile <- file.path(outpath,'seqtab_noprior.csv');
write.csv(seqtabfinal,outfile);

outfile <- file.path(outpath, 'seqtab_noprior.rds')
saveRDS(seqtabfinal, file=outfile)

ASV_sequences <- colnames(seqtabfinal);
print("number of ASVs")
print(length(ASV_sequences))
sprintf("%s - before writing out ASVs", Sys.time())

fasta_file= file.path(outpath,'ASV_sequences.fasta');
library(seqinr);
write.fasta(sequences=as.list(ASV_sequences), names=sprintf("ASV_temp_%d",1:length(ASV_sequences)), fasta_file)

seqtabfinal_counts_raw = as.data.frame(seqtabfinal) #Name table column as `ASV_temp_#` instead of full ASV sequence.
colnames(seqtabfinal_counts_raw) = sprintf("ASV_temp_%d",1:ncol(seqtabfinal_counts_raw))
library(reshape2);
sprintf("%s - before melting",Sys.time())
seqtabfinal_counts_raw$oligos_id = rownames(seqtabfinal);
seqtabfinal_counts = melt(seqtabfinal_counts_raw, id.vars="oligos_id", variable.name="asv_temp_id",value.name="count" )
seqtabfinal_counts = seqtabfinal_counts[seqtabfinal_counts$count!=0,];
sprintf("%s - before writing out asv_counts.csv",Sys.time())

count_table_outfile <- file.path(outpath, "asv_counts.csv");
write.csv(seqtabfinal_counts, count_table_outfile, row.names = F);

sprintf("%s - before writing out error profiles",Sys.time())
outfile <- file.path(outpath, 'errR.rds')
saveRDS(errR, file=outfile)

outfile <- file.path(outpath, 'errF.rds')
saveRDS(errF, file=outfile)

outfile <- file.path(outpath, 'derepFs.rds')
saveRDS(derepFs, file=outfile)

outfile <- file.path(outpath, 'derepRs.rds')
saveRDS(derepRs, file=outfile)

sprintf("%s - ending script",Sys.time())
