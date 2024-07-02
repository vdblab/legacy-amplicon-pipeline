import glob
import os
envvars:
    "HOME"


configfile: os.path.join(workflow.basedir, "config/config.yaml")


def get_files_from_input_dir(config):
    # we dont print by default without the verbose option  cause this issue:
    # https://github.com/snakemake/snakemake/issues/1297
    # print("getting input fastqs from " + config["input_directory"])
    files = glob.glob(config["input_directory"] + "/*f*q*")
    if not [x for x in files if "R1_001.f" in x]:
        raise ValueError("Pattern 'R1_001.f' not found in input files; please double check and/or fix names")
    # print(files)
    oligos_files = glob.glob(config["input_directory"] + "/*.oligos")
    #print(oligos_files)
    files_dict = {
        "readsf": [x for x in files if "R1_001.f" in x],
        "readsr": [x for x in files if "R2_001.f" in x],
        "oligos": oligos_files[0]
    }
    for k,v in files_dict.items():
        #print("   :-:", k, ": ", v)
        if v == None:
            raise ValueError("Unable to find required {} file in input directory".format(k))
    return files_dict

input_files = get_files_from_input_dir(config)



def get_files_from_db_dir(config):
    #print("getting paths for database files in " + config["db_directory"])
    db_base = "bacteria_and_archea.16SrRNA.fna"
    anno_base = "bacteria_and_archea.16SrRNA.id_and_taxonomy_v2.txt"
    files_dict = {
        "db":   glob.glob(config["db_directory"] + "/*" + db_base)[0],
        "anno": glob.glob(config["db_directory"] + "/*" + anno_base)[0]
    }
    for k,v in files_dict.items():
        #print("   :-:", k, ": ", v)
        if v == None:
            raise ValueError("Unable to find required {} file in db directory".format(k))
    return files_dict

db_files = get_files_from_db_dir(config)


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
            samples.append(id_pool.replace(thisid, cleanid))


onstart:
    print("Samples found in oligos file:\n- " + "\n- ".join(samples))

rule all:
    input:
        "dada2_results/blast_out/ASV_annotated.txt",
        "dada2_results/ASV_sequences.fasta",
        "dada2_results/asv_counts.csv",
        "dada2_results/blast_out/blast_passed_not_passed.txt"

rule merge_lanes_if_needed:
    input:
        readsf = input_files["readsf"],
        readsr = input_files["readsr"]
    output:
        temp("reads1.fastq.gz"),
        temp("reads2.fastq.gz")
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
        readsf="reads1.fastq",
        readsr="reads2.fastq",
        barcodes="barcodes.fastq",
    container: "docker://ghcr.io/vdblab/biopython:1.70a"
    log: "logs/primer_removal.log"
    message: "01 - removing primer sequences from fastq pools"
    params:
        primerf=config['primerf'],
        primerr=config['primerr']
    script: "scripts/strip_addons3_py3.py"


rule clean_oligos:
    input: input_dir = config["input_directory"]
    output: oligos=os.path.basename(input_files["oligos"]) + ".clean"
    container: "docker://ghcr.io/vdblab/dada2legacy:1.18.0"
    message: "02 - cleaning oligos file"
    log: "logs/oligo_cleaning.log"
    script: "scripts/oligos_cleaning_cl.R"


rule convert_oligos_to_mapping_file:
    input:
        oligosfile=os.path.basename(input_files["oligos"]) + ".clean"
    output:
        "1.map.txt",
        "2.map.txt"
    container: "docker://ghcr.io/vdblab/biopython:1.70a"
    message: "03 - converting oligos file to mapping file for demultiplexing"
    log: "logs/convert_oligos.log"
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
    log: "logs/guess_encoding.log"
    script: "scripts/guest-encoding.py"

rule add_demultiplex_info_to_fastq:
    input:
        reads1="reads1.fastq",
        reads2="reads2.fastq",
        map1="1.map.txt",
        map2="2.map.txt",
        encoding="encoding.txt",
        barcodes="barcodes.fastq"
    output:
        labelled_reads_F = "demultiplex_r1/seqs.fastq",
        labelled_reads_R = "demultiplex_r2/seqs.fastq"
    params:
        demultiplex_r1= "demultiplex_r1",
        demultiplex_r2= "demultiplex_r2"
    log: "logs/add_demultiplex_info.log"
    container: "docker://ghcr.io/vdblab/qiime:1.9.1"
#    conda: "envs/py2.yaml"
    message: "05 - labelling reads with sample id"
    shell: """
    echo "spliting out forward libraries" >> {log}
    split_libraries_fastq.py -i {input.reads1} -o {params.demultiplex_r1} -b {input.barcodes} --store_demultiplexed_fastq -q 0 -s 100000000 -m {input.map1} --barcode_type 12 --max_barcode_errors 4 -n 5000 -p 0.00001 -r 1000000 --phred_offset `cut {input.encoding} -f2`
    echo "splitting out reverse libraries" >> {log}
    split_libraries_fastq.py -i {input.reads2} -o {params.demultiplex_r2} -b {input.barcodes} --store_demultiplexed_fastq -q 0 -s 100000000 -m {input.map2} --barcode_type 12 --max_barcode_errors 4 -n 5000 -p 0.00001 -r 1000000 --phred_offset `cut {input.encoding} -f2`

   rm reads1.fastq
   rm reads2.fastq
   rm scrap.fastq
   echo "done" >> {log}
   """

rule set_up_dummy_files:
    input:
        map1="1.map.txt",
    params:
        sampleid_dir= "sampleids"
    output:
        sampleid_files= expand("sampleids/{i}.sample", i=samples),
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
    container: "docker://ghcr.io/vdblab/qiime:1.9.1"
#    conda: "envs/py2.yaml"
    message: "06 - extracting reads for sample"
    shell: """
    sampleid=$(cat {input.sampleid_file})
    temp_file=temp"$sampleid".txt
    echo $sampleid > $temp_file

    if [ ! -d isolated_oligos ]; then mkdir isolated_oligos; fi
    a=`date`
    echo "$a: Extracting sampleid $sampleid" >> {log}
    filter_fasta.py -f {input.reads1} --sample_id_fp {input.sampleid_file} -o isolated_oligos/${{sampleid}}_R1.fastq 2>> {log}
    filter_fasta.py -f {input.reads2} --sample_id_fp {input.sampleid_file} -o isolated_oligos/${{sampleid}}_R2.fastq 2>> {log}
    gzip -f isolated_oligos/${{sampleid}}_R1.fastq 2>> {log}
    gzip -f isolated_oligos/${{sampleid}}_R2.fastq 2>> {log}
    rm $temp_file
    touch tmp
    """


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
    threads: 8 # this is hardcoded in dadascript_noprior_ag.R.  For now...
    message: "07 - Identifying ASVs with DADA2"
    #container: "docker.pkg.github.com://vdblab/vdblab-pipelines/16s-r:548aae2c1be691d2dde659fe38e4efe0d7849906"
    # shell: """
    # ls
    # R -e 'packageVersion("dada2"); packageVersion("ShortRead")'
    # Rscript AG_16S/Docker/dadascript_noprior_ag.R isolated_oligos
    # """
    script: "scripts/dadascript_noprior_ag.R"


rule run_ASV_annotation:
    input:
        fasta = "dada2_results/ASV_sequences.fasta",
        db = db_files["db"],
        annotation = db_files["anno"]
    output:
        annotated_asvs="dada2_results/blast_out/ASV_annotated.txt",
        blast = "dada2_results/blast_out/blast_passed_not_passed.txt",
    message: "08 - Annotating ASVs with BLAST"
    log: "logs/blast_asvs.log"
    container: "docker://ghcr.io/vdblab/legacyampliconblast:0.0.1"
    script: "scripts/blast_simplified_dada2_cl.R"
