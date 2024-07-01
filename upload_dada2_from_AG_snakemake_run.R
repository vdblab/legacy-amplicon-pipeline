#Oct/30/2019
#
#It uploads DADA2 dada to database.
#It takes as input the folder name that contains fastq.gz files
#It uses the structure:
#     isolated_oligos/ and isolated_oligos/blast_out to get corresponding files
#File names:
#fasta_file: (...)/isolated_oligos/ASV_sequences.fasta
#counts_file: (...)/isolated_oligos/counts_asv.csv
#ASV_annotated: isolated_oligos/blast_out/ASV_annotated.txt
#blast_details: isolated_oligos/blast_out/blast_passed_not_passed.txt
#
#Sep/23/2019
#
#I changed this script to take as input the folders that contains pool.fasta files.
#Then it identifies otu_call_closed/ and the corresponding files.
#
#
#Oct242018
#I changed the script to take a file with list of alpha_diversity or otu_table as input.
#This scripts uploads alpha_diversity and counts data to database from command line!
#It is used inside the cluster to get it automatically piped after running 16S pipeline.

# 2022-03-10 NRW simplify, use vdbR, remove absolute paths, remove extraneous code, remove trailing semicolons,
# match with snakemake output structure

library(seqinr)
library(stringr)
library(biomformat)
library(vdbR)
connect_database()

args = commandArgs(trailingOnly = T)



get_data_from_query_OTU = function(query_id, ...) {

  query = '';

  args = list(...);

  if(query_id==-1){
    #query is input in args[1];

    query = args[[1]];
  }

  if(query_id==0) {

    table_name = args[1];
    q = postgresqlReadTable(psql_con, name=table_name); #Get `table` directly, without need of a query.

  }

  if(query_id==0.1){
    #get column names from a table;
    table_name = args[1];

    query= sprintf("select column_name, data_type  from information_schema.columns where table_name='%s'",table_name);
  }

  if(query_id==0.2){
    #get maximum `key` value from a table;
    table_name = args[1];

    query= sprintf("select max(key) from %s",table_name);
  }

  if(query_id==0.3){
    #get counts from a table;
    table_name = args[1];

    query= sprintf("select count(*) from %s",table_name);
  }

  if(query_id==0.5 | query_id=="showtables"){
    query="SELECT tablename FROM pg_catalog.pg_tables where tablename like '%_ag' order by tablename";
  }

  if(query_id==1) {
    #Get all distinct `(sample_id, source) and total counts` from count table.

    query = "select cag.sampleid, cag.source_table, sum(cag.count) as total_coverage, count(cag.otu_key) as n_otus from counts_ag as cag group by cag.sampleid, cag.source_table";

    #table_name = "samples";
    #q = get_data_from_query_OTU(0,table_name);

  }
  if(query!=''){
    q = dbGetQuery(psql_con, query);
  }
}

upload_data_from_query_OTU_check_and_submission <- function(table_name, d_set_to_upload ){
  #Check if d_set is in proper format to be uploaded in `table_name`.

  #Add uploaded_date to date_frame.
  uploaded_date = format(Sys.time(),"%m-%d-%Y");
  d_set_to_upload$uploaded_date = uploaded_date;

  #Assign incremental key values starting on `maximum` value in current table.

  q_key_max_cur = get_data_from_query_OTU(0.2,table_name);
  if(is.na(q_key_max_cur$max)){
    q_key_max_cur$max = 0;
  }
  d_set_to_upload$key = ( 1:length(d_set_to_upload$uploaded_date) ) + q_key_max_cur$max;

  #Get columns for `column names`
  q_column_names = get_data_from_query_OTU(0.1,table_name);

  #stop("Add step to compare columns, re-order data, add key and upload!");
  if( ! all(is.element(q_column_names$column_name, colnames(d_set_to_upload) )) ){
    stop( sprintf("Some fields in table %s are not defined in upload dataset", table_name));
  }

  d_set_to_upload_reordered = d_set_to_upload[,q_column_names$column_name]

  #Add the step to write the data into table!
  print("data getting uploaded!");

  dbWriteTable(con, table_name, value = d_set_to_upload_reordered, append = TRUE, row.names = FALSE);
  #my_dbWriteTable(con, table_name, d_set_to_upload_reordered);
  print( sprintf("upload: done! %d rows", length(d_set_to_upload_reordered$key)));

}


update_data_from_query_OTU_check_and_submission <- function(table_name, d_set_to_upload){
  #This version check duplicates in table (according to table keys) and add only new, unique data.
  #It intermediates this step with a tempory table.
  #Some interesting discussion on skip duplicates: https://stackoverflow.com/questions/1009584/how-to-emulate-insert-ignore-and-on-duplicate-key-update-sql-merge-with-po
  #Rationale:
  #     1. Create temp_updating table with same structure as target table
  #     2. Upload data to temporary table
  #     3. Retrieve data in temporary table that is not in target table
  #     4. Remove temporary table.
  #     5. Return only novel data to be uploaded;


  print("");
  print(sprintf("table: %s", table_name));
  #Clean temp_table
  temp_table="temp_updating";
  query_check_temp = sprintf("SELECT EXISTS (SELECT 1 FROM   information_schema.tables WHERE table_name = '%s')",
                             temp_table);
  q_check = dbGetQuery(con,query_check_temp);
  if(q_check$exists){
    dbSendQuery(con, sprintf("drop table %s", temp_table));
  }

  #Add uploaded_date to date_frame.
  uploaded_date = format(Sys.time(),"%m-%d-%Y");
  d_set_to_upload$uploaded_date = uploaded_date;

  #Assign incremental key values starting on `maximum` value in current table.
  q_key_max_cur = get_data_from_query_OTU(0.2,table_name);
  if(is.na(q_key_max_cur$max)){
    q_key_max_cur$max = 0;
  }
  d_set_to_upload$key = ( 1:length(d_set_to_upload$uploaded_date) ) + q_key_max_cur$max;

  #Get columns for `column names`
  q_column_names = get_data_from_query_OTU(0.1,table_name);
  #q_column_unique_names = get_data_from_query_OTU(0.1,sprintf("%s_unique",table_name));

  #stop("Add step to compare columns, re-order data, add key and upload!");
  if( ! all(is.element(q_column_names$column_name, colnames(d_set_to_upload) )) ){
    stop( sprintf("Some fields in table %s are not defined in upload dataset", table_name));
  }

  d_set_to_upload_reordered = d_set_to_upload[,q_column_names$column_name]

  print("    creating temporary table to check for duplicates");
  query_create_temp_table = sprintf("CREATE TABLE %s AS SELECT * FROM %s WHERE 1=2",temp_table, table_name);
  dbSendQuery(con,query_create_temp_table);

  dbWriteTable(con, temp_table, value = d_set_to_upload_reordered, append = TRUE, row.names = FALSE);
  #my_dbWriteTable(con, table_name, d_set_to_upload_reordered);

  ###getting data that does not violate unique constrain in target table.
  query_select_unique_key = paste("select ccu.column_name, ccu.constraint_name from information_schema.table_constraints as tc",
                                  " natural join information_schema.constraint_column_usage as ccu ",
                                  #" on tc.table_name=ccu.table_name and tc.constraint_name=ccu.constraint_name ",
                                  sprintf(" where tc.table_name='%s' and tc.constraint_type='UNIQUE'",table_name),
                                  sep = "");
  q_unique_key = dbGetQuery(con, query_select_unique_key);
  if(length(unique(q_unique_key$constraint_name))>1){
    stop("This table seems to have more than one rule for unique constraint. This script is not ready to work with it");
  }
  unique_constraint_str = paste(q_unique_key$column_name,collapse = ",");
#  query_select_novel_data = sprintf("select * from %s where %s as not in (select %s from %s)",
#                                    temp_table,
#                                    unique_constraint_str,
#                                    unique_constraint_str,
#                                    table_name);
  unique_constraint_join_str = paste(sprintf("t1.%s=t2.%s",q_unique_key$column_name,q_unique_key$column_name), collapse = " AND ")

  query_select_novel_data = paste(sprintf("select t1.* from %s as t1 left join %s as t2", temp_table, table_name),
                                  sprintf(" on %s ", unique_constraint_join_str),
                                  sprintf("where t2.%s is NULL", q_unique_key$column_name[1]),
                                  collapse = " ");

  d_set_unique = dbGetQuery(con,query_select_novel_data);

  print(sprintf("%d/%d novel rows to be updated!",
                dim(d_set_unique)[1],
                dim(d_set_to_upload_reordered)[1]) );
  if(length(d_set_unique)==0){
    print("      No updates! No novel data.");
    return();
  }
  upload_data_from_query_OTU_check_and_submission(table_name,d_set_unique);
}

upload_data_from_query_OTU <- function(query_number, ...){

  args = list(...)
  input_data_file = args[[1]]


  if(query_number==33.1){

    table_name = "asv_sequences_ag";

    asv_fasta_file = input_data_file;
    asv_fasta = read.fasta(asv_fasta_file,as.string = T);

    if( anyDuplicated(tolower(unlist(getSequence(asv_fasta,as.string = T))))){
      stop("I did not expected all ASVs to be unique.");
    }

    d_set=data.frame(asv_sequence= tolower(unlist(getSequence(asv_fasta,as.string = T))));
          #1
    update_data_from_query_OTU_check_and_submission(table_name, d_set);
  }

  if(query_number==33.2){
    table_name = "asv_counts_ag";

    asv_fasta_file = args[[1]]
    asv_counts_file = args[[2]]

    asv_fasta = read.fasta(asv_fasta_file, as.string = T);
    dt_fasta = data.table(asv_sequence = factor(tolower(unlist(getSequence(asv_fasta,as.string = T)))),
                          asv_temp_id=names(asv_fasta));

    q0_asv_sequences = get_data_from_query_OTU(0, "asv_sequences_ag");
    q0_asv_sequences$asv_sequence = factor(q0_asv_sequences$asv_sequence);
    if(!all(dt_fasta$asv_sequence %in% q0_asv_sequences$asv_sequence)){
      stop("I expected all asv_sequence to to in database.")
    }

    dt_fasta_with_key = merge(dt_fasta,
                              q0_asv_sequences,
                              by="asv_sequence");

    asv_counts = fread(asv_counts_file)

    asv_counts_and_asv_key = merge(dt_fasta_with_key[,.(asv_temp_id,asv_key=key)],
                                   asv_counts[,.(asv_temp_id=asv_temp_id,
                                                 oligos_id,
                                                 sampleid=sapply(oligos_id, function(x) {strsplit(x,"\\.\\.")[[1]][1]}),
                                                 count)]);

    d_set=d_set=data.frame(asv_key=asv_counts_and_asv_key$asv_key,
                           sampleid=asv_counts_and_asv_key$sampleid,
                           oligos_id=asv_counts_and_asv_key$oligos_id,
                           count=asv_counts_and_asv_key$count);


    d_set=d_set=data.frame(asv_key=asv_counts_and_asv_key$asv_key,
                           sampleid=asv_counts_and_asv_key$sampleid,
                           oligos_id=asv_counts_and_asv_key$oligos_id,
                           count=asv_counts_and_asv_key$count);

    #Check whether oligoid is already in database.
    q0_asv_counts = get_data_from_query_OTU(0,"asv_counts_ag");
    if(length(q0_asv_counts)==0){
      ind_sampleid_in_database = rep(F,length(d_set$oligos_id));
    }else{
      oligos_id_already_in_database=unique( q0_asv_counts[, .(oligos_id)][[1]])
      #u_sampleid_and_run_date_to_be_uploaded=unique( d_set[, c("sampleid", "run_date")] );
      ind_sampleid_in_database = d_set$oligos_id %in% oligos_id_already_in_database;
    }


    if(all(ind_sampleid_in_database)){
      print("All oligo_ids are already uploaded to asv_counts_ag");
      return();
    }
    if(any(ind_sampleid_in_database)){
      warning("some sampleids are already in database. Uploading only new ones");
      d_set = d_set[!ind_sampleid_in_database,];
    }

    update_data_from_query_OTU_check_and_submission(table_name, d_set);

  }

 if(query_number==33.21){
    table_name = "asv_alpha_diversity_ag";

    asv_fasta_file = args[[1]]
    asv_counts_file = args[[2]]

    asv_fasta = read.fasta(asv_fasta_file, as.string = T);
    dt_fasta = data.table(asv_sequence = factor(tolower(unlist(getSequence(asv_fasta,as.string = T)))),
                          asv_temp_id=names(asv_fasta));

    q0_asv_sequences = get_data_from_query_OTU(0, "asv_sequences_ag");
    q0_asv_sequences$asv_sequence = factor(q0_asv_sequences$asv_sequence);
    if(!all(dt_fasta$asv_sequence %in% q0_asv_sequences$asv_sequence)){
      stop("I expected all asv_sequence to to in database.")
    }

    dt_fasta_with_key = merge(dt_fasta,
                              q0_asv_sequences,
                              by="asv_sequence");

    asv_counts = fread(asv_counts_file)

    asv_counts_and_asv_key = merge(dt_fasta_with_key[,.(asv_temp_id,asv_key=key)],
                                   asv_counts[,.(asv_temp_id=asv_temp_id,
                                                 oligos_id,
                                                 sampleid=sapply(oligos_id, function(x) {strsplit(x,"\\.\\.")[[1]][1]}),
                                                 count)]);

    d_set_raw=d_set=data.table(asv_key=asv_counts_and_asv_key$asv_key,
                           sampleid=asv_counts_and_asv_key$sampleid,
                           oligos_id=asv_counts_and_asv_key$oligos_id,
                           count=asv_counts_and_asv_key$count);

    simpson_r_f <- function(x) {
      p_i = x/sum(x);
      r = return( 1/(sum(p_i^2)) );
    }

    shannon_f <- function(x) {
      p_i = x/sum(x);
      r = return( sum(-p_i*log(p_i)) );
    }

    d_set = d_set_raw[,.(simpson_reciprocal=simpson_r_f(count),
                         shannon= shannon_f(count), count_total=sum(count)),
                      by=.(sampleid,oligos_id)];

    path_pool = dirname(asv_fasta_file);

    d_set$path_pool = path_pool;

    update_data_from_query_OTU_check_and_submission(table_name, as.data.frame(d_set));
  }


if(query_number==33.3){
    table_name = "asv_annotation_blast_ag";

    asv_fasta_file = args[[1]]
    annotation_file = args[[2]]

    asv_fasta = read.fasta(asv_fasta_file, as.string = T);
    dt_fasta = data.table(asv_sequence = factor(tolower(unlist(getSequence(asv_fasta,as.string = T)))),
                          asv_temp_id=names(asv_fasta));
    q0_asv_sequences = get_data_from_query_OTU(0, "asv_sequences_ag");
    q0_asv_sequences$asv_sequence = factor(q0_asv_sequences$asv_sequence);
    if(!all(dt_fasta$asv_sequence %in% q0_asv_sequences$asv_sequence)){
      stop("I expected all asv_sequence to to in database.")
    }
    dt_fasta_with_key = merge(dt_fasta,
                              q0_asv_sequences,
                              by="asv_sequence");

    dt_annotation = fread(annotation_file);

    dt_annotation_and_key = merge(dt_fasta_with_key[,.(asv_temp_id,asv_key=key)],
                                  dt_annotation,
				  by="asv_temp_id");
    dt_annotation_and_key=dt_annotation_and_key[order(asv_key)]


    split_annotation = function(annotation_str_raw){
      annotation_str = sapply(annotation_str_raw, function(x){strsplit(x,";")[[1]]});
      df_annot = as.data.frame(sapply(annotation_str, function(x) {strsplit(x,split = "__")}));
      field_names = t(df_annot)[,1];
      field_values = t(df_annot)[,2];
      names(field_values) = field_names
      return(field_values)
    }

    y = lapply(as.character(dt_annotation_and_key$annotation),
                           function(x) split_annotation(x) );

    df_annotation_parsed = (Reduce("rbind", y));

    d_set = data.frame(asv_key = dt_annotation_and_key$asv_key,
                       kingdom = df_annotation_parsed[,"k"],
                       phylum  = df_annotation_parsed[,"p"],
                       class  = df_annotation_parsed[,"c"],
                       ordr  = df_annotation_parsed[,"o"],
                       family = df_annotation_parsed[,"f"],
                       genus = df_annotation_parsed[,"g"],
                       species = df_annotation_parsed[,"s"],
                       blast_pass = dt_annotation_and_key$pident); #97% identify,

    update_data_from_query_OTU_check_and_submission(table_name, d_set);


  }

  if(query_number==33.31){
    table_name = "asv_annotation_blast_detailed_ag";

    asv_fasta_file = args[[1]]
    annotation_details_file = args[[2]]

    asv_fasta = read.fasta(asv_fasta_file, as.string = T);
    dt_fasta = data.table(asv_sequence = factor(tolower(unlist(getSequence(asv_fasta,as.string = T)))),
                          asv_temp_id=names(asv_fasta));
    q0_asv_sequences = get_data_from_query_OTU(0, "asv_sequences_ag");
    q0_asv_sequences$asv_sequence = factor(q0_asv_sequences$asv_sequence);
    if(!all(dt_fasta$asv_sequence %in% q0_asv_sequences$asv_sequence)){
      stop("I expected all asv_sequence to to in database.")
    }
    dt_fasta_with_key = merge(dt_fasta,
                              q0_asv_sequences,
                              by="asv_sequence");

    dt_annotation_details = fread(annotation_details_file);

    dt_annotation_and_key = merge(dt_fasta_with_key[,.(asv_temp_id,asv_key=key)],
                                  dt_annotation_details,
                                  by.x="asv_temp_id",
                                  by.y="ASVId");
    dt_annotation_and_key=dt_annotation_and_key[order(asv_key)]

    d_set_input = as.data.frame(dt_annotation_and_key)

    d_set=data.frame(asv_key=d_set_input$asv_key,
		     accession=d_set_input$accession,
                     name=d_set_input$name,
                     unique_name=d_set_input$unique_name,
                     query_length=d_set_input$query_length,
                     align_length=d_set_input$align_length,
                     pident=d_set_input$pident,
                     nident=d_set_input$nident,
                     score=d_set_input$score,
                     length_ratio=d_set_input$length_ratio,
                     pident_97=d_set_input$pident_97,
                     length_ratio_1_1=d_set_input$length_ratio_1.1,
                     passed=d_set_input$passed)

    update_data_from_query_OTU_check_and_submission(table_name, d_set);


  }


}


input_folder_cur = args[[1]]

print(sprintf("input_folder_cur: %s",input_folder_cur));



file_set = c(asv_fasta_file = file.path(input_folder_cur, "dada2_results/ASV_sequences.fasta"),
             asv_counts_file = file.path(input_folder_cur, "dada2_results/asv_counts.csv"),
             asv_annotated_file = file.path(input_folder_cur, "dada2_results/blast_out/ASV_annotated.txt"),
             asv_blast_details_file = file.path(input_folder_cur, "dada2_results/blast_out/blast_passed_not_passed.txt"))

if(! all(file.exists(file_set) ) ){
    stop(sprintf("path %s does not have all DADA2 files", input_folder_cur))
}

print(sprintf("uploading the results from %s", input_folder_cur))

upload_data_from_query_OTU(33.1,
                           file_set["asv_fasta_file"])
upload_data_from_query_OTU(33.2,
                           file_set["asv_fasta_file"],
                           file_set["asv_counts_file"])
upload_data_from_query_OTU(33.21,
                           file_set["asv_fasta_file"],
                           file_set["asv_counts_file"])
upload_data_from_query_OTU(33.3,
                           file_set["asv_fasta_file"],
                           file_set["asv_annotated_file"])
upload_data_from_query_OTU(33.31,
                           file_set["asv_fasta_file"],
                           file_set["asv_blast_details_file"])
