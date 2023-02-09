from collections import defaultdict
from pathlib import Path

import pandas as pd


rule index_bam:
    input:
        bam=INPUT_DIR+"/{sample}.bam"
    output:
        bam_index=INPUT_DIR+"/{sample}.bam.bai"
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -c 'pysam.index({input.bam})'"


rule count_junctions:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=rules.index_bam.output.bam_index,
        known="known_SJ"
    output:
        junctions=OUTPUT_DIR+"/J1/{sample}.J1.gz",
        library_stats=OUTPUT_DIR + "/J1/{sample}.library_stats.txt"
    threads: THREADS
    params:
        primary=("", "-p")[config["primary"]],
        unique=("", "-u")[config["unique"]]
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.count_junctions "
        "-i {input.bam} "
        "-k {input.known} "
        "-o {output.junctions} "
        "-l {output.library_stats} "
        "{params.primary} {params.unique} "
        "-t {threads}"


checkpoint gather_library_stats:
    input:
        library_stats=expand("{out}/J1/{sample}.library_stats.txt", out=OUTPUT_DIR, sample=samples)
    output:
        tsv=OUTPUT_DIR+"/aggregated_library_stats.tsv"
    shell:
        "python3 -m workflow.scripts.gather_library_stats "
        "{OUTPUT_DIR}/J1  "
        "-o {output.tsv}"


rule aggregate_junctions:
    input:
        junctions=rules.count_junctions.output.junctions,
        library_stats=rules.count_junctions.output.library_stats
    output:
        aggregated_junctions=OUTPUT_DIR+"/J2/{sample}.J2.gz"
    params:
        min_offset=config["min_offset"],
        min_intron_length=config["min_intron_length"],
        max_intron_length=config["max_intron_length"]
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.aggregate_junctions "
        "-i {input.junctions} "
        "-s {input.library_stats} "
        "-o {output.aggregated_junctions} "
        "--min_offset {params.min_offset} "
        "--min_intron_length {params.min_intron_length} "
        "--max_intron_length {params.max_intron_length}"


def get_org(s):
    stats = pd.read_table(checkpoints.gather_library_stats.get().output.tsv, index_col=0)
    return stats.loc[s, "genome"]


rule annotate_junctions:
    input:
        aggregated_junctions=rules.aggregate_junctions.output.aggregated_junctions,
        genome=lambda wildcards: f"genomes/{get_org(wildcards.sample)}.fa",
        known_sj=lambda wildcards: f"known_SJ/{get_org(wildcards.sample)}.ss.tsv.gz"
    output:
        annotated_junctions=OUTPUT_DIR+"/J3/{sample}.J3.gz"
    conda: "../envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.annotate_junctions "
         "-i {input.aggregated_junctions} "
         "-k {input.known_sj} "
         "-f {input.genome} "
         "-o {output.annotated_junctions}"


rule choose_strand:
    input:
        annotated_junctions=rules.annotate_junctions.output.annotated_junctions,
        ranked_list=lambda wildcards: f"known_SJ/{get_org(wildcards.sample)}.ranked.txt"
    output:
        junction_stats=OUTPUT_DIR+"/J4/{sample}.junction_stats.txt",
        stranded_junctions=OUTPUT_DIR+"/J4/{sample}.J4.gz"
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.choose_strand "
        "-i {input.annotated_junctions} "
        "-r {input.ranked_list} "
        "-o {output.stranded_junctions} "
        "-s {output.junction_stats}"


checkpoint gather_junction_stats:
    input:
        junction_stats=expand("{out}/J4/{sample}.junction_stats.txt", out=OUTPUT_DIR, sample=samples)
    output:
        tsv=OUTPUT_DIR+"/aggregated_junction_stats.tsv"
    run:
        d = defaultdict(list)
        for replicate in input.junction_stats:
            p = Path(replicate)
            name = Path(p.stem).stem
            with p.open("r") as f:
                d["replicate"].append(name)
                for line in f:
                    if line.startswith("-"):
                        break
                    left, right = line.strip().split(": ")
                    d[left].append(right)
        df = pd.DataFrame(d)
	if not df.empty:
        	df.sort_values(by="replicate")
        df.to_csv(output.tsv, index=False, sep="\t")


rule filter_junctions:
    input:
        stranded_junctions=rules.choose_strand.output.stranded_junctions
    output:
        filtered_junctions=OUTPUT_DIR+"/J6/{sample}.J6.gz"
    params:
        entropy=config["entropy"],
        total_count=config["total_count"],
        gtag=("", "-g")[config["gtag"]]
    conda: "../envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.filter "
         "-i {input.stranded_junctions} "
         "-e {params.entropy} "
         "-c {params.total_count} "
         "{params.gtag} "
         "-o {output.filtered_junctions}"


rule merge_junctions:
    input:
         stranded_junctions=expand("{out}/J4/{sample}.J4.gz", sample=samples, out=OUTPUT_DIR)
    output:
         merged_junctions=OUTPUT_DIR+"/J4/merged_junctions.J4.gz"
    conda: "../envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.merge_junctions "
         "{input.stranded_junctions} "
         "-o {output.merged_junctions}"