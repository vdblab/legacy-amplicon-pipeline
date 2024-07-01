#!/bin/bash
#
#Script to run dada2 pipeline and blast annotation
#
#dadascript_noprior_ag.R
#	-Dada2 script will output:
#		- asv_sequences.fasta
#		- asv_counts.csv
#			-in long format.
#
#Blast_simplified_dada2_cl:
#
#	- Annotate asv_sequences according to blast 16S database;
#
#TO DO:
#	- Annotate using GTDB
#
#	- DATABASE table:
#		- asv_counts
#			- sampleid, oligoid
#		- asv_sequences
#		- asv_annotation_blast
#		- asv_annotation_gtdb
#
#		(...)
#

set -e;

path_with_fastq_file=$1;


path_demultiplexed="$path_with_fastq_file""/isolated_oligos/";

echo "Running pipeline on this path: $path_demultiplexed"

Rscript dadascript_noprior_ag.R $path_demultiplexed;

fasta_file="$path_demultiplexed""ASV_sequences.fasta";

echo "starting blast" >&2

Rscript blast_simplified_dada2_cl.R $fasta_file
