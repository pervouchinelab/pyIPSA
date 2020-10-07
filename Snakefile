import pysam

configfile: "config/config.yaml"
INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]
THREADS = config["threads"]
samples, = glob_wildcards(INPUT_DIR + "/{sample}.bam")


rule all:
    input:
         expand("{out}/{sample}/1.j.gz", sample=samples, out=OUTPUT_DIR),
         expand("{out}/{sample}/stats.txt", sample=samples, out=OUTPUT_DIR)


rule count_junctions:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=INPUT_DIR+"/{sample}.bam.bai"
    output:
        junctions=OUTPUT_DIR+"/{sample}/1.j.gz",
        stats=OUTPUT_DIR+"/{sample}/stats.txt"
    params:
        threads = THREADS
    shell:
        "python -m workflow.scripts.junctions_counter "
        "-i {input.bam} "
        "-s known_splice_junctions "
        "-o {output.junctions} "
        "-a {output.stats} "
        "-t {params.threads}"


rule index_bam:
    input:
        bam=INPUT_DIR+"/{sample}.bam"
    output:
        bam_index=INPUT_DIR+"/{sample}.bam.bai"
    run:
        pysam.index(input.bam)
