import pysam
import requests
import os
configfile: "config.yaml"
samples, = glob_wildcards(config["INPUT_DIR"]+"/{sample}.bam")

if not samples:
    r = requests.get("https://cb.skoltech.ru/dp/pyipsa/miniPCAWG.bam", allow_redirects=True)
    if not os.path.exists(config["INPUT_DIR"]):
        os.makedirs(config["INPUT_DIR"])
    open(config["INPUT_DIR"] + "/miniPCAWG.bam", "wb").write(r.content)
    samples = ["miniPCAWG"]

rule all:
    input:
        expand("{output_dir}/{sample}/annotated.ssj.gz", sample=samples, output_dir=config["OUTPUT_DIR"]),
        expand("{output_dir}/{sample}/unprocessed.ssc.gz", sample=samples, output_dir=config["OUTPUT_DIR"])

rule annotate_splice_site_junctions:
    input:
        aggssj=config["OUTPUT_DIR"]+"/{sample}/aggregated.ssj.gz",
        genome="genomes/hg19.fa",
        introns="introns/hg19.introns"
    output:
        annssj=config["OUTPUT_DIR"]+"/{sample}/annotated.ssj.gz"
    shell:
        "python annotate.py -ssj {input.aggssj} -o {output.annssj} -fa {input.genome} -a {input.introns}"


rule aggregate_splice_site_junctions:
    input:
        ssj=config["OUTPUT_DIR"]+"/{sample}/unprocessed.ssj.gz",
        log=config["OUTPUT_DIR"]+"/{sample}/library_analysis.txt",
    output:
        aggssj=config["OUTPUT_DIR"]+"/{sample}/aggregated.ssj.gz"
    params:
        min_offset=config["MIN_OFFSET"]
    shell:
        "python aggregate.py -i {input.ssj} -o {output.aggssj} -l {input.log} -m {params.min_offset}"

rule analyze_library:
    input:
        bam=config["INPUT_DIR"]+"/{sample}.bam",
        bam_index=config["INPUT_DIR"]+"/{sample}.bam.bai",
        ssj=config["OUTPUT_DIR"]+"/{sample}/unprocessed.ssj.gz"
    output:
        log=config["OUTPUT_DIR"]+"/{sample}/library_analysis.txt"
    shell:
        "python library_analyzer.py -i {input.bam} -ssj {input.ssj} > {output.log}"

rule count_splice_site_continuous:
    input:
        bam=config["INPUT_DIR"]+"/{sample}.bam",
        bam_index=config["INPUT_DIR"]+"/{sample}.bam.bai",
        ssj=config["OUTPUT_DIR"]+"/{sample}/unprocessed.ssj.gz"
    output:
        ssc=config["OUTPUT_DIR"]+"/{sample}/unprocessed.ssc.gz"
    shell:
        "python ssc_counter.py -i {input.bam} -ssj {input.ssj} | sort -k1,1 -k2,2n | gzip > {output.ssc}"

rule count_splice_site_junctions:
    input:
        bam=config["INPUT_DIR"]+"/{sample}.bam",
        bam_index=config["INPUT_DIR"]+"/{sample}.bam.bai"
    output:
        ssj=config["OUTPUT_DIR"]+"/{sample}/unprocessed.ssj.gz"
    shell:
        "python ssj_counter.py -i {input.bam} | sort -k1,1 -k2,2n | gzip > {output.ssj}"

rule index_bam:
    input:
        bam=config["INPUT_DIR"]+"/{sample}.bam"
    output:
        bam_index=config["INPUT_DIR"]+"/{sample}.bam.bai"
    run:
        pysam.index(input.bam)

rule download_annotation:
    output:
        gtf="annotations/hg19.gtf"
    shell:
        """
        wget -O {output.gtf}.gz {config[URL_ANNOTATIONS][hg19]}
        gunzip {output.gtf}.gz
        """

rule download_genome:
    output:
        genome="genomes/hg19.fa"
    shell:
        """
        wget -O {output.genome}.gz {config[URL_GENOMES][hg19]}
        gunzip {output.genome}.gz
        """