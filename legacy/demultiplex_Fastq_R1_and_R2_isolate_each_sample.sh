#May/22/2019
#This version outputs an individual file for each sampleid in oligos file.
#
#It uses qiime command: `XXXXXX` and output files to `(...)/isolated_oligos/`
#
#Script to create fastq demultiplexed file
#May/18/2018
#I modified this script to demultiplex R1 and R2 separately!
#Jan/24/2018.

#Activate qiime1 environment

#echo "DON'T USE THIS SCRIPT, demultiplexing R1 and R2 files separetely generates files of different length. If you really want to use it, check how to use it!”
#exit;

#source ~/.source_miniconda
# export PATH=$PATH:/home/gomesa/miniconda3/bin/
set -e

data_path=$1;
#list_of_sampleid_to_keep=$2;

data_path=`echo "$data_path" | awk '{if(substr($1,length($1),1)!="/") {$1=$1"/";} print $1}'`;

cur_path=`pwd`;
echo "your starting path is $cur_path";

log_file="temp_log"; #_"$time_stamp;

echo "running folder: $data_path" >> $log_file;

Rscript oligos_cleaning_cl.R $data_path;

n_R1_files=`ls $data_path*R1*fastq.gz | wc -l`;
n_R2_files=`ls $data_path*R2*fastq.gz | wc -l`;
n_oligos_files=`ls $data_path*.oligos | wc -l`;

echo "n_R1_files : $n_R1_files";
echo "n_R2_files : $n_R2_files";
if [[ $n_R1_files -ne 1 ]] || [[ $n_R2_files -ne 1 ]] || [[ $n_oligos_files -ne 1 ]]
then
echo "n_R1_files : $n_R1_files";
echo "n_R2_files : $n_R2_files";
output_message=$'program terminated: Unmatch in number of files: \n\t n_R1_files=$n_R1_files \n\t n_R2_files=$n_R2_files \n\t n_oligos_files=$n_oligos_files';
echo $output_message;
echo $output_message >> $log_file;
exit;
fi;

cd "$data_path";

#source activate qiime1;

R1_fastq_gz_files=`ls *R1*fastq.gz`;
R2_fastq_gz_files=`ls *R2*fastq.gz`;
oligos_file=`ls *.oligos`;

file_match_test=`echo "$R1_fastq_gz_files $R2_fastq_gz_files" | awk '{R1_file=$1; R2_file=$2; gsub("_R1_","_R2_",R1_file); if(R1_file!=R2_file) {out="Error:file mismatch\n" "R1_file="R1_file "\nR2_file="R2_file "\n"} else {out=1}; print out }'`;

#echo "check 1"

echo "check $file_match_test";

if [[ ! $file_match_test == 1 ]]
then
	echo "Error:Mismatch in R1 and R2 filenames”;
	echo "Error:Mismatch in R1 and R2 filenames”  >> $log_file;
	exit;
fi;

echo "check 1.5"

#Run demultiplex using strip_addons.py. I believe this command creates a barcodes.fastq file, necessary for `split_libraries_fastq.py` and reads1.fastq as well as reads2.fastq files. I believe `extract_barcodes.py` from qiime assumes barcode is at beginning of read and strip_addons fix this assumption.

#a="strip_addons.py $R1_fastq_gz_files $R2_fastq_gz_files -remove_bar_primer -pdiffs 1 -fw_primer AYTGGGYDTAAAGNG -rev_primer CCGTCAATTYHTTTRAGT"

a="python strip_addons3.py $R1_fastq_gz_files $R2_fastq_gz_files -remove_bar_primer -fw_primer AYTGGGYDTAAAGNG -rev_primer CCGTCAATTYHTTTRAGT -pdiffs 1"

echo $a;

#echo "I deactivated strip_addons.py part for some checks. Make sure it is reactivated"
$a;

echo "check 2"

convert_oligos_to_qiime_mapping_file_rj.py $oligos_file;

guest_encoding="/home/gomesa/projects/human_microbiota/library_cluster/guest-encoding.py"
phred_value=`head -1000 reads1.fastq | awk 'NR%4==0' | $guest_encoding | awk '$0~/phred33/ {a=33} $0~/phred64/ {a=a""64} {print a}' `;

#split_libraries_fastq.py -i reads1.fastq -o test_r1/ -b barcodes.fastq --store_demultiplexed_fastq -q 0 -s 100000000 -m 1.map.txt --barcode_type 24 --max_barcode_errors 4
split_libraries_fastq.py -i reads1.fastq -o demultiplex_r1/ -b barcodes.fastq --store_demultiplexed_fastq -q 0 -s 100000000 -m 1.map.txt --barcode_type 12 --max_barcode_errors 4 -n 5000 -p 0.00001 -r 1000000 --phred_offset $phred_value;

split_libraries_fastq.py -i reads2.fastq -o demultiplex_r2/ -b barcodes.fastq --store_demultiplexed_fastq -q 0 -s 100000000 -m 2.map.txt --barcode_type 12 --max_barcode_errors 4 -n 5000 -p 0.00001 -r 1000000 --phred_offset $phred_value;

rm reads1.fastq
rm reads2.fastq
rm scrap.fastq


#List of sampleids
awk 'NR>1 {print $1}' 1.map.txt > list_of_sampleids.txt

#Filtering each sampleid to its own file:
#Dec/04/2019
#Following this instruction to make parallel demultiplex: https://stackoverflow.com/questions/17307800/how-to-run-given-function-in-bash-in-parallel/17316302
demultiplex_each_sample() {
	sampleid=$1;
	temp_file=temp"$sampleid".txt
	echo $sampleid > $temp_file;
	a=`date`;
	echo "$a: Extracting sampleid $sampleid" >> temp_isolated_log.txt;
	output_file_r1="isolated_oligos/""$sampleid""_R1.fastq";
	output_file_r2="isolated_oligos/""$sampleid""_R2.fastq";
	filter_fasta.py -f demultiplex_r1/seqs.fastq --sample_id_fp $temp_file -o $output_file_r1;
	filter_fasta.py -f demultiplex_r2/seqs.fastq --sample_id_fp $temp_file -o $output_file_r2;
	gzip -f $output_file_r1;
	gzip -f $output_file_r2;
	rm $temp_file;

}
export -f demultiplex_each_sample;

if [ ! -d isolated_oligos ]; then mkdir isolated_oligos; fi
#mkdir isolated_oligos;
a=`date`
echo "$a: Starting sampleid extraction";
echo "$a: Starting sampleid extraction" > temp_isolated_log.txt;

parallel -a list_of_sampleids.txt demultiplex_each_sample

if [ 1 -eq 0 ]; then
echo "CODE NOT RUN TO TEST PARALLEL FUNCTION"
for sampleid in `cat list_of_sampleids.txt`; do
echo $sampleid > temp.txt;
date >> temp_isolated_log.txt
echo "Starting sampleid $sampleid" >> temp_isolated_log.txt;
output_file_r1="isolated_oligos/""$sampleid""_R1.fastq";
output_file_r2="isolated_oligos/""$sampleid""_R2.fastq";
filter_fasta.py -f demultiplex_r1/seqs.fastq --sample_id_fp temp.txt -o $output_file_r1;
filter_fasta.py -f demultiplex_r2/seqs.fastq --sample_id_fp temp.txt -o $output_file_r2;
gzip -f $output_file_r1;
gzip -f $output_file_r2;
done;
fi

rm -rf demultiplex_r1/
rm -rf demultiplex_r2/

source deactivate;

cd "$cur_path";
