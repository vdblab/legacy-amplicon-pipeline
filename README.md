## `AG_16S` pipeline

### Config Parameters

- `input_directory`: path to directory with input files
- `db_directory`: path to directory containing NCBI's 16S database. Currently at `/data/brinkvd/resources/`
- `trunclen`: desired length of the reads after trimming.  the original AG pipeline was hardcoded to 180bp, but for data generated outside of our center this should be changed accordingly.
- `minasvlen`: minimum length of ASV; maximum is set to 2x `trunclen`.  This is used to filter ASVs.  The original pipeline defulted to 300, with a max of 360 (2*180)
- `ncores`:  number of cores to give to the DADA2 rule


### Usage
```
snakemake --cores 4 --use-conda --use-singularity --snakefile ./AG_16S/Snakefile --config input_directory=/mnt/disk1/test_data/ db_directory=/mnt/disk1/16S_blast/ --directory output
```

For example, submitting a job via snakemake would be:

```
snakemake --jobs 32 --use-conda  --use-singularity --snakefile AG_16S/Snakefile --config input_directory=/home/daia1/my_workdir/other_pipeline/dada2/test/ trunclen=180 minasvlen=300 ncores=32  db_directory=/data/brinkvd/resources/ --directory results/ --cluster 'bsub -n 32 -R "rusage[mem=4]" -W 12:00 -e pipeline.err -o pipeline.out'
```
