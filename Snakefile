import glob
import gzip
import os
from contextlib import redirect_stderr
import traceback
import pandas as pd



configfile: os.path.join(workflow.basedir, "config/config.yaml")


# skip a schema for now:
for req in ["R1", "R2", "oligos", "blast_db", "blast_anno"]:
    assert req in config, f"config argument {req} is mandatory"

input_files =  {
        "readsf": config["R1"],
        "readsr": config["R2"],
        "oligos": config["oligos"]
    }


db_files = {
        "db":   config["blast_db"],
        "anno": config["blast_anno"]
    }



samp_id_taboo_chars = '[_%+;: -]+'

samples = []
with open (input_files["oligos"], "r") as inf:
    for line in inf:
        if line.startswith("BARCODE"):
            id_pool = line.split()[3]
            if ".." not in id_pool:
                print(f"Line with misspecified oligo_id:\n{line}")
                raise ValueError("To maintain compatibility with previous pipeline runs, the oligos file must be specified as {sampleid}..{pool} in the 4th column.  Exiting")
            thisid = id_pool.split("..")[0]
            # borrowed this regex from the barcode cleaning script
            cleanid = re.sub(samp_id_taboo_chars,'.',thisid)
            thissample = id_pool.replace(thisid, cleanid)
            assert thissample not in samples, "Error: duplicate samples detected in oligos file"
            samples.append(thissample)

# we add "Unassigned" to samples in some outputs
# see http://qiime.org/scripts/split_libraries_fastq.html

localrules:
   all,
   output_manifest,

onstart:
    print("Samples found in oligos file:\n- " + "\n- ".join(samples))

rule all:
    input:
        "dada2_results/blast_out/ASV_annotated.txt",
        "dada2_results/ASV_sequences.fasta",
        "dada2_results/asv_counts.csv",
        "dada2_results/blast_out/blast_passed_not_passed.txt",
        "multiqc/multiqc_report.html",
        "demux/manifest.tsv",
	"isolated_oligos/Unassigned_R1.fastq.gz",
	"isolated_oligos/Unassigned_R2.fastq.gz",

rule merge_lanes_if_needed:
    input:
        readsf = input_files["readsf"],
        readsr = input_files["readsr"]
    output:
        temp("reads1.fastq.gz"),
        temp("reads2.fastq.gz")
    threads: 1
    shell:
        """
        cat {input.readsf} >> reads1.fastq.gz
        cat {input.readsr} >> reads2.fastq.gz
        """

rule remove_primers:
    input:
        readsf = "reads1.fastq.gz",
        readsr = "reads2.fastq.gz"
    output:
        readsf=temp("reads1.fastq"),
        readsr=temp("reads2.fastq"),
        barcodes="barcodes.fastq",
    container: "docker://ghcr.io/vdblab/biopython:1.70a"
    log: "logs/primer_removal.log"
    message: "01 - removing primer sequences from fastq pools"
    resources:
        mem_mb=lambda wc, attempt: 12 * 1024 * attempt,
        runtime=lambda wc, attempt: 4 * 60 * attempt,
    threads: 1
    params:
        primerf=config['primerf'],
        primerr=config['primerr']
    script: "scripts/strip_addons3_py3.py"

# rule clean_oligos:
#     input: oligos = config["oligos"]
#     output: oligos=os.path.basename(input_files["oligos"]) + ".clean"
#     container: "docker://ghcr.io/vdblab/dada2legacy:1.18.0"
#     message: "02 - cleaning oligos file"
#     threads: 1
#     log: "logs/oligo_cleaning.log"
#     script: "scripts/oligos_cleaning_cl.R"


rule convert_oligos_to_mapping_file:
    input:
        oligosfile=input_files["oligos"]
    output:
        "1.map.txt",
        "2.map.txt"
    container: "docker://ghcr.io/vdblab/biopython:1.70a"
    message: "03 - converting oligos file to mapping file for demultiplexing"
    log: "logs/convert_oligos.log"
    threads: 1
    params:
        outdir=".",
        primerf=config['primerf'],
        primerr=config['primerr']
    script: "scripts/convert_oligos.py"


rule guess_encoding_of_fastq:
    input: "reads1.fastq"
    container: "docker://ghcr.io/vdblab/biopython:1.70a"
    output: "encoding.txt"
    message: "04 - determining the encoding of the FASTQ quality scores"
    threads: 1
    log: "logs/guess_encoding.log"
    script: "scripts/guess-encoding.py"

rule add_demultiplex_info_to_fastq:
    input:
        reads="reads1.fastq",
        map="1.map.txt",
        encoding="encoding.txt",
        barcodes="barcodes.fastq"
    output:
        labelled_reads = "demultiplex_r1/seqs.fastq",
	fna = temp("demultiplex_r1/seqs.fna"),
    params:
        demultiplex= "demultiplex_r1",
    container: "docker://ghcr.io/vdblab/qiime:1.9.1"
    resources:
        mem_mb=lambda wc, attempt: 12 * 1024 * attempt,
        runtime=lambda wc, attempt: 4 * 60 * attempt,
    message: "05 - labelling reads with sample id"
    shell: """
    split_libraries_fastq.py -i {input.reads} -o {params.demultiplex} -b {input.barcodes} --store_demultiplexed_fastq -q 0 -s 100000000 -m {input.map} --barcode_type 12 --max_barcode_errors 4 -n 5000 -p 0.00001 -r 1000000 --retain_unassigned_reads --phred_offset `cut {input.encoding} -f2`
    """

use rule add_demultiplex_info_to_fastq as add_demultiplex_info_to_fastq_R2 with:
    input:
        reads="reads2.fastq",
        map="2.map.txt",
        encoding="encoding.txt",
        barcodes="barcodes.fastq"
    output:
        labelled_reads = "demultiplex_r2/seqs.fastq",
	fna = temp("demultiplex_r2/seqs.fna"),
    params:
        demultiplex= "demultiplex_r2"

rule set_up_dummy_files:
    input:
        map1="1.map.txt",
    params:
        sampleid_dir= "sampleids"
    output:
        sampleid_files= expand("sampleids/{i}.sample", i=samples),
	unassigned_file = "sampleids/Unassigned.sample"
    message: "05.5 - splitting up map file for distributed demultiplexing"
    log: "logs/split_map.log"
    script: "scripts/split_up_map_file.py"

rule demultiplex:
    input:
        reads1="demultiplex_r1/seqs.fastq",
        reads2="demultiplex_r2/seqs.fastq",
        sampleid_file="sampleids/{sampleid}.sample"
    output:
        reads1="isolated_oligos/{sampleid}_R1.fastq.gz",
        reads2="isolated_oligos/{sampleid}_R2.fastq.gz"

    log: "logs/demultiplex_{sampleid}.log"
    threads: 1
    container: "docker://ghcr.io/vdblab/qiime:1.9.1"
    message: "06 - extracting reads for sample"
    resources:
        mem_mb=16 * 1024,
        runtime = lambda wc, attempt: 4 * 60 * attempt,
    shell: """
    sampleid=$(cat {input.sampleid_file})
    temp_file=temp"$sampleid".txt
    echo $sampleid > $temp_file

    a=`date`
    echo "$a: Extracting sampleid $sampleid" >> {log}
    filter_fasta.py -f {input.reads1} --sample_id_fp {input.sampleid_file} -o isolated_oligos/${{sampleid}}_R1.fastq 2>> {log}
    filter_fasta.py -f {input.reads2} --sample_id_fp {input.sampleid_file} -o isolated_oligos/${{sampleid}}_R2.fastq 2>> {log}
    gzip -f isolated_oligos/${{sampleid}}_R1.fastq 2>> {log}
    gzip -f isolated_oligos/${{sampleid}}_R2.fastq 2>> {log}
    rm $temp_file
    touch tmp
    """

use rule demultiplex as demultiplex_unassigned with:
    input:
        reads1="demultiplex_r1/seqs.fastq",
        reads2="demultiplex_r2/seqs.fastq",
        sampleid_file="sampleids/Unassigned.sample"
    output:
        reads1="isolated_oligos/Unassigned_R1.fastq.gz",
        reads2="isolated_oligos/Unassigned_R2.fastq.gz"

    log: "logs/demultiplex_Unassigned.log"

# https://stackoverflow.com/questions/37874936/how-to-check-empty-gzip-file-in-python
def gz_size(fname):
    with gzip.open(fname, "rb") as f:
        return f.seek(0, whence=2)

def write_manifest_and_missing(
    sample_ids, fastq_template, manifest_path, missing_path, paired=True
):
    R1 = []
    R2 = []
    for s in sample_ids:
        R1_path = fastq_template.format(sample=s, dir=1)
        if os.path.exists(R1_path) and gz_size(R1_path) > 0:
            R1.append(R1_path)
        else:
            R1.append("")

        R2_path = fastq_template.format(sample=s, dir=2)
        if os.path.exists(R2_path) and gz_size(R2_path) > 0:
            R2.append(R2_path)
        else:
            R2.append("")
    manifest = pd.DataFrame({"sample_id": sample_ids, "R1": R1, "R2": R2})
    unassigned_manifest = pd.DataFrame({
        "sample_id": "Unassigned",
        "R1":[fastq_template.format(sample="Unassigned", dir=1)],
	"R2": [fastq_template.format(sample="Unassigned", dir=2)],
	})
    manifest = manifest[["sample_id", "R1", "R2"]]
    unassigned_manifest = unassigned_manifest[["sample_id", "R1", "R2"]]
    if paired:
        is_incomplete = (manifest["R1"] == "") | (manifest["R2"] == "")
    else:
        is_incomplete = manifest["R1"] == ""

    pd.concat([manifest[is_incomplete], unassigned_manifest]).to_csv(missing_path, sep="\t", index=False)

    manifest = manifest[~is_incomplete]
    manifest.to_csv(manifest_path, sep="\t", index=False)

rule output_manifest:
    """
    Output a file with sample ids extracted from the oligo files,
    and demultiplexed R1 and R1 fastq files
    """
    input:
        expand(
            "isolated_oligos/{sample}_R{dir}.fastq.gz",
            sample=samples,
            dir=[1, 2],
        ),
	"isolated_oligos/Unassigned_R1.fastq.gz",
    output:
        manifest=f"demux/manifest.tsv",
        missing=f"demux/missing_or_incomplete.tsv",
    log:
        e=f"log/output_manifest.e",
    run:
        with open(log.e, "w") as ef, redirect_stderr(ef):
            try:
                fastq_template = os.path.join(
                    os.getcwd(), "isolated_oligos/{sample}_R{dir}.fastq.gz"
                )
                write_manifest_and_missing(
                    samples, fastq_template, output["manifest"], output["missing"]
                )
            except Exception as e:
                traceback.print_exc(file=ef)


rule run_dada2:
    # the input_fastqs acts as a gather signal, waiting for all the fastqs to be available
    input:
#        input_fastqs = expand("isolated_oligos/{i}_R2.fastq.gz", i=samples),
        readsf=expand("isolated_oligos/{i}_R1.fastq.gz", i=samples),
        readsr=expand("isolated_oligos/{i}_R2.fastq.gz", i=samples),
    params:
        input_dir = "isolated_oligos",
        trunclen = config["trunclen"],
        outdir_name = "dada2_results",
        minasvlen= config["minasvlen"],
        ncores = 16
    output:
        fasta = "dada2_results/ASV_sequences.fasta",
        counts = "dada2_results/asv_counts.csv"
    container: "docker://ghcr.io/vdblab/dada2legacy:1.18.0"
    log: "logs/dada2script.log"
    threads: 16
    resources:
        mem_mb=48 * 1024,
        runtime = lambda wc, attempt: 12 * 60 * attempt,
    message: "07 - Identifying ASVs with DADA2"
    script: "scripts/dadascript_noprior_ag.R"

rule sample_fastqc_report:
    input:
        readsf="isolated_oligos/{sample}_R1.fastq.gz",
        readsr="isolated_oligos/{sample}_R2.fastq.gz",
    output:
        report_R1="fastqc_reports/{sample}_R1_fastqc.html",
        report_R2="fastqc_reports/{sample}_R2_fastqc.html",
        zip_R1="fastqc_reports/{sample}_R1_fastqc.zip",
        zip_R2="fastqc_reports/{sample}_R2_fastqc.zip",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://staphb/fastqc:0.11.9"
    threads: 2
    resources:
        mem_mb=4 * 1024,
    shell:
        """
        fastqc \
            --outdir {params.outdir} \
            --threads {threads} \
            --noextract \
            {input}
        """

rule multiqc:
    input:
        readsf=expand("fastqc_reports/{sample}_R1_fastqc.html", sample=samples),
    output:
        "multiqc/multiqc_report.html"
    container: "docker://ewels/multiqc:v1.12"
    shell:"""
            multiqc        -o multiqc \
            --force --filename multiqc_report.html ./
    """



rule run_ASV_annotation:
    input:
        fasta = "dada2_results/ASV_sequences.fasta",
        db = db_files["db"],
        annotation = db_files["anno"]
    output:
        annotated_asvs="dada2_results/blast_out/ASV_annotated.txt",
        blast = "dada2_results/blast_out/blast_passed_not_passed.txt",
    message: "08 - Annotating ASVs with BLAST"
    threads: 8 # hardcoded for now
    resources:
        mem_mb=16 * 1024,
        runtime = lambda wc, attempt: 4 * 60 * attempt,
    log: "logs/blast_asvs.log"
    container: "docker://ghcr.io/vdblab/legacyampliconblast:0.0.1"
    script: "scripts/blast_simplified_dada2_cl.R"
