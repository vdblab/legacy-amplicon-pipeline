#Dec/25/2019
#
#Clean oligos file:
#
#1. trim `.txt` from oligos files.
#2. Change "\r" to "\n"
#3. Add "..pool###" to castori center oligoid.
#4. Crash if there is duplicated oligoid.
#5. Creates a bkp file (.pre_cleaning.bkp) only if it does not pre-exist.
#

#data_path="~/projects/human_microbiota/data/16S/dada2/temp_Sample_pooltest_complete/"

#args = commandArgs(trailingOnly = T);

#data_path=args[[1]];


log <- file(snakemake@log[[1]], open="wt")
sink(log)
data_path = snakemake@input[["input_dir"]]


if(data_path=="./"){
  #It needs the full path name to extract the "pool###" from the path.
  data_path =getwd()
}
oligos_file=dir(data_path, pattern = "\\.oligos$|\\.oligos\\.txt$", full.names=TRUE)

print(oligos_file)

#print(sprintf("N oligos: %d", n_oligos_files))
if(length(oligos_file) != 1){
  error_str = sprintf("Path %s should have exactly 1 oligos file.",data_path);
  #print(error_str)
  stop(error_str)
}


oligos_file_no_txt = outdir_path = snakemake@output[["oligos"]]

#if(oligos_file!=oligos_file_no_txt){
##   print("removing .txt ending for oligos file.");
##   file.rename(oligos_file,oligos_file_no_txt);
## }else{
##   print("File ends in .oligos");
## }

#awk '$0' was necessary to clean up white space that occured in some files.
system(sprintf("tr '\r' '\n' < %s | awk '$0' > temp_%s; mv temp_%s %s",
               oligos_file,
               basename(oligos_file_no_txt),
               basename(oligos_file_no_txt),
               oligos_file_no_txt));

dt_oligos_primer_lines = read.table(oligos_file,header = F,nrows = 2,
                                    sep = "\t",comment.char = "#",
                                    fill = T);

if(is.null(dt_oligos_primer_lines$V4)){
  dt_oligos_primer_lines$V4="";
}

dt_oligos_primer_lines$V4[is.na(dt_oligos_primer_lines$V4)]="";

dt_oligos_primer_lines[is.na(dt_oligos_primer_lines)]="";

print("primer lines done!")
dt_oligos = read.table(oligos_file,header = F,skip = 2,sep = "\t",comment.char = "#")
dt_oligos[is.na(dt_oligos)]=""

if(grepl("Sample_",data_path) & grepl("_complete",data_path)){
  #Rename pools from Castori center. Folder format: Sample_POOL_complete
  pool_info = strsplit(data_path,"Sample_")[[1]][2];
  pool_info = strsplit(pool_info,"_complete")[[1]][1];

  #Only add "..pool###" when it is not already in row and it is not a comment.
  #Some files have more info about it and I am keeping them.
  ind_missing_pool_and_not_comment = (!grepl(pattern = pool_info, as.character(dt_oligos$V4)) & !grepl("\\#",dt_oligos$V1));

  if(any(grepl("pool", as.character(dt_oligos$V4)[ind_missing_pool_and_not_comment]))){
    stop("Oligoid seem to belong to a different `pool###` (path and name mismatch)")
  }
  if(any(ind_missing_pool_and_not_comment) ){
    dt_oligos$V4=as.character(dt_oligos$V4);
    dt_oligos$V4[ind_missing_pool_and_not_comment] = sapply(as.character(dt_oligos$V4[ind_missing_pool_and_not_comment]),
                          function(x) {strsplit(x,"\\.\\.")[[1]][1]});
    dt_oligos$V4[ind_missing_pool_and_not_comment] = paste(dt_oligos$V4[ind_missing_pool_and_not_comment],sprintf("..%s",pool_info),sep="");
  }else{
    print("All oligoids in `.oligos` file already contained proper `..pool###`");
  }

}

if(any(duplicated(dt_oligos$V4))){
   stop("Duplicated oligoids were not expected!");
}
oligos_file_bkp = sprintf("%s.pre_cleaning.bkp",oligos_file_no_txt);

## if(file.exists(oligos_file_bkp)){
##   #bkp_str = "You have already performed cleaning. Make sure to remove .oligos.bkp file to unlock and run it again"
##   print("backup string already exists"); #oligos_file_bkp);
## #  stop(bkp_str);
## }else{
## 	file.rename(oligos_file_no_txt,oligos_file_bkp);
## }

write.table(dt_oligos_primer_lines,oligos_file_no_txt,row.names = F,col.names = F,sep="\t",quote = F);
write.table(dt_oligos,oligos_file_no_txt,row.names = F,col.names = F,sep="\t",append = T,quote = F)
