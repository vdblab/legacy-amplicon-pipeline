## Legacy pipeline


### Environment Prerequisites

- Snakemake
- Singularity
- A snakemake profile configured and the `SNAKEMAKE_PROFILE` variable set

### Config Parameters

- `R1`: a list of paths to the R1 reads (eg `[/path/to/pool_L001_R1.fq.gz,/path/to/pool_L002_R1.fq.gz]`)
- `R2`: a list of paths to the R2 reads (eg `[/path/to/pool_L001_R2.fq.gz,/path/to/pool_L002_R2.fq.gz]`)
- `oligos`: path to oligos file
- `trunclen`: desired length of the reads after trimming.  the original AG pipeline was hardcoded to 180bp, but for data generated outside of our center this should be changed accordingly.
- `minasvlen`: minimum length of ASV; maximum is set to 2x `trunclen`.  This is used to filter ASVs.  The original pipeline defulted to 300, with a max of 360 (2*180)
- `blast_db`: path to the fasta in the blast db folder
- `blast_anno`: path to the fasta in the blast db taxonomy annotation file

### Usage
```
snakemake  --config trunclen=180 minasvlen=300 R1=[$PWD/test/test_input/test_R1_001.fastq.gz] R2=[$PWD/test/test_input/test_R2_001.fastq.gz] oligos=$PWD/test/test_input/pool1059.oligos  blast_db=/data/brinkvd/resources_old/bacteria_and_archea.16SrRNA.fna blast_anno=/data/brinkvd/resources_old/bacteria_and_archea.16SrRNA.id_and_taxonomy_v2.txt --directory $PWD/output
```
