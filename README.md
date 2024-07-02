## `AG_16S` pipeline

### Config Parameters

- `input_directory`: path to directory with input files
- `db_directory`: path to directory containing NCBI's 16S database. Currently at `/data/brinkvd/resources/`
- `trunclen`: desired length of the reads after trimming.  the original AG pipeline was hardcoded to 180bp, but for data generated outside of our center this should be changed accordingly.
- `minasvlen`: minimum length of ASV; maximum is set to 2x `trunclen`.  This is used to filter ASVs.  The original pipeline defulted to 300, with a max of 360 (2*180)
- `ncores`:  number of cores to give to the DADA2 rule


### Usage
```
snakemake --cores 4 --use-conda --use-singularity  --config trunclen=180 minasvlen=300 input_directory=$PWD/test/test_input/ db_directory=/data/brinkvd/resources_old/ --directory $PWD/output
```
