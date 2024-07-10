library(data.table)

if (exists("snakemake")){
                                        #print(str(snakemake))
    print("parsing snakemake args")
    # see https://stackoverflow.com/questions/64101921/
    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    fasta_file <- snakemake@input[["fasta"]]
    output_path <- sprintf("%s/blast_out/", dirname(fasta_file));#args[2];
    database_file = snakemake@input[["db"]]
    annotation_file = snakemake@input[["annotation"]]

} else{
    stop("Script can only be run via Snakemake")
}

args = snakemake@input

re_run_boolean=T

dir.create(output_path)

print("Parsing annotation file")
dt_annotation <- fread(annotation_file, sep="\t", header=T)
dt_annotation <- as.data.frame(dt_annotation)

#The annotation file will map tax_id_species to a taxonomic annotation.
if(! all(colnames(dt_annotation) == c("accession_id", "tax_id_species", "taxon"))){
  stop("I expected dt_annotation to have the columns: accession_id, tax_id_species, taxon");
}

rownames(dt_annotation) = dt_annotation$accession_id;



library(plyr)
library(ape)

blast.function <- function(fasta_file, output_path, re_run=T){
  #fasta<-x
  #dir.fasta<-paste(miseq.dir,fasta,sep="")
  print(output_path);
#  setwd(output_path);
  out_results = paste(output_path,"results.txt",sep="");

  if(re_run){

    if (file.exists(out_results)) file.remove(out_results)

    #provide the directory info for blastn unix executable
    cmds<-paste("blastn -db", database_file," -max_target_seqs 500  -num_threads 8 -query", fasta_file, "-outfmt \"6 qseqid staxids saccver stitle qlen length nident pident bitscore evalue score\" -out", out_results, sep=" ")
    #print(cmds);
    #Sep/26/2018: I added `max_target_seqs` to make sure we query enough sequences. There is risk of bias using a small number of target sequences: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166
    sapply(cmds, system)  # takes ~30 min
  }
  #if error, may need to manually delete results.txt

  results <- fread(out_results, sep="\t",header=FALSE,stringsAsFactors=FALSE)
  setnames(results, c("ASVId", "taxid", "accession", "species", "query_length", "align_length", "nident", "pident", "bitscore", "evalue", "score"))

  results$tax_id_species = dt_annotation[as.character(results$accession),"tax_id_species"];
  #ASVs<-data.frame(ASVId=unique(results$ASVId),stringsAsFactors=FALSE)
  ASVs<-unique(results$ASVId)

  ##following identifies the accessions with the best match, including tie scores
  blast_result_extract <- function(x){
    #x<-ASVs[1]
    #print(x)
    results_ASV<-results[results$ASVId==x,]
    #find the tie score hits
    results_ASV_top<-results_ASV[results_ASV$bitscore==results_ASV$bitscore[1],]

    #for the species, replace the first space with underscore, delete second space and beyond, then alphabetical order
    results_ASV_top$species_short<-sub(" ","_",results_ASV_top$species)
    results_ASV_top$species_short<-sub(" .*$","",results_ASV_top$species_short)
    results_ASV_top<-results_ASV_top[order(results_ASV_top$species_short),]

    #identify unique names, put them in alphabetical order, and merge them, and also add the % identity
    name<-paste(paste(unique(results_ASV_top$species_short),collapse=";"),results_ASV$align_length[1],results_ASV$pident[1],sep=";")

    #generate a unique name with the oldest taxid
    results_ASV_top<-results_ASV_top[order(as.numeric(results_ASV_top$tax_id_species)),]
    old_taxid<-results_ASV_top$taxid[1]
    accession<-results_ASV_top$accession[1];
    unique_name<-results_ASV_top$species_short[1]
    align_length<-results_ASV$align_length[1]
    pident<-results_ASV$pident[1]
    score<-results_ASV$score[1]
    evalue<-results_ASV$evalue[1];
    query_length<-results_ASV$query_length[1]
    nident<-results_ASV$nident[1]
    blast_results<-data.frame(ASVId=x,taxid=old_taxid, accession=accession, name=name,unique_name=unique_name,query_length=query_length, align_length=align_length,pident=pident,nident=nident,score=score,evalue=evalue)
    return(blast_results)

  }

  blasted<-ldply(ASVs, function(x) blast_result_extract(x));

  return(blasted)
  #setwd(miseq.dir)
}


print("executing blast")
blasted <- blast.function(fasta_file = fasta_file, output_path=output_path,re_run=re_run_boolean)
print("Filtering blast results")
#1
#this next section selects ASVs that should be added back if removed during chimera selection
#see how much of the ASV sequence BLAST was able to align
blasted$length_ratio<-blasted$query_length/blasted$align_length

#examine ASVs with a pident score 97% or better
blasted$pident_97<-blasted$pident >=97
blasted_97<-blasted[blasted$pident_97==TRUE,]
#plot the results
#hist(blasted_97$length_ratio)
#it's a bimodal distribution, so let's use a cut-off of 1.1

blasted$length_ratio_1.1<-blasted$length_ratio <= 1.1
blasted$passed<-blasted$length_ratio_1.1 == TRUE & blasted$pident_97 == TRUE

#so these are the ASVs that blasted well enough to avoid chimera checking
blasted_passed<-blasted[blasted$passed==TRUE,];
blasted_not_passed<-blasted[blasted$passed==FALSE,];

#if(grep(".fasta",fasta_file)){
#	output_file_base = gsub(".fasta","",fasta_file);
#}else{
#	output_file_base = fasta_file;
#}
output_file_base = output_path; #paste(output_path, gsub(".fasta","",basename(fasta_file)),sep="");

save_blast_passed = sprintf("%sblast_passed.txt",output_file_base);
save_blast_not_passed = sprintf("%sblast_not_passed.txt", output_file_base);

save_blast_detailed = sprintf("%sblast_passed_not_passed.txt",output_file_base);

#The output of blast is send to this file to be load in database.
write.table(blasted_passed, save_blast_passed, sep="\t", row.names=F);
write.table(blasted_not_passed, save_blast_not_passed, sep="\t", row.names=F);
write.table(rbind(blasted_passed,blasted_not_passed), save_blast_detailed, sep="\t", row.names=F);


if(length(blasted_passed$ASVId)){
  ASV_annotated_p1 <- data.frame( asv_temp_id = blasted_passed$ASVId,
                                  blast_pass=T,
                                  accession=blasted_passed$accession,
                                  annotation=dt_annotation[as.character(blasted_passed$accession),"taxon"],
                                  tax_id_species=dt_annotation[as.character(blasted_passed$accession),"tax_id_species"],
                                  pident=blasted_passed$pident,
                                  evalue=blasted_passed$evalue,
                                  score=blasted_passed$score)
}
if(length(blasted_not_passed$ASVId)){
  ASV_annotated <- data.frame( asv_temp_id = blasted_not_passed$ASVId,
                               blast_pass=F,
                               accession=blasted_not_passed$accession,
                               annotation=dt_annotation[as.character(blasted_not_passed$accession),"taxon"],
                               tax_id_species=dt_annotation[as.character(blasted_not_passed$accession),"tax_id_species"],
                               pident=blasted_not_passed$pident,
                               evalue=blasted_not_passed$evalue,
                               score=blasted_not_passed$score);
if(length(blasted_passed$ASVId)){
ASV_annotated = rbind(ASV_annotated_p1, ASV_annotated)
}
}else{
ASV_annotated = ASV_annotated_p1;
}

ASV_annotated_file  = sprintf("%sASV_annotated.txt", output_file_base);

write.table(ASV_annotated, ASV_annotated_file,sep="\t", row.names=F);
