rule index_bam:
    input:
        bam=INPUT_DIR+"/{sample}.bam"
    output:
        bam_index=INPUT_DIR+"/{sample}.bam.bai"
    run:
        pysam.index(input.bam)


checkpoint count_junctions:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=rules.index_bam.output.bam_index,
        splice_sites="known_SJ"
    output:
        junctions=OUTPUT_DIR+"/J1/{sample}.J1.gz",
        stats=OUTPUT_DIR+"/J1/{sample}.stats"
    params:
        threads=THREADS,
        primary=("", "-p")[config["primary"]],
        unique=("", "-u")[config["unique"]]
    shell:
        "python3 -m workflow.scripts.count_junctions "
        "-i {input.bam} "
        "-s {input.splice_sites} "
        "-o {output.junctions} "
        "-l {output.stats} "
        "{params.primary} {params.unique} "
        "-t {params.threads}"


rule aggregate_junctions:
    input:
        junctions=OUTPUT_DIR+"/J1/{sample}.J1.gz",
        stats=OUTPUT_DIR+"/J1/{sample}.stats"
    output:
        aggregated_junctions=OUTPUT_DIR+"/J2/{sample}.J2.gz"
    params:
        min_offset=config["min_offset"],
        strand_mode=config["strand_mode"],
        min_intron_length=config["min_intron_length"],
        max_intron_length=config["max_intron_length"]
    shell:
        "python3 -m workflow.scripts.aggregate_junctions "
        "-i {input.junctions} "
        "-s {input.stats} "
        "-o {output.aggregated_junctions} "
        "--strand_mode {params.strand_mode} "
        "--min_offset {params.min_offset} "
        "--min_intron_length {params.min_intron_length} "
        "--max_intron_length {params.max_intron_length}"


def get_org(s):
    with open(checkpoints.count_junctions.get(sample=s).output.stats) as file:
        return file.readline().strip().split(" is ")[1]


rule annotate_junctions:
    input:
        aggregated_junctions=rules.aggregate_junctions.output.aggregated_junctions,
        genome=lambda wildcards: f"genomes/{get_org(wildcards.sample)}.fa",
        known_sj=lambda wildcards: f"known_SJ/{get_org(wildcards.sample)}.ss.tsv.gz"
    output:
        annotated_junctions=OUTPUT_DIR+"/J3/{sample}.J3.gz"
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
        stranded_junctions=OUTPUT_DIR+"/J4/{sample}.J4.gz"
    shell:
        "python3 -m workflow.scripts.choose_strand "
        "-i {input.annotated_junctions} "
        "-r {input.ranked_list} "
        "-o {output.stranded_junctions}"


rule filter_junctions:
    input:
        stranded_junctions=rules.choose_strand.output.stranded_junctions
    output:
        filtered_junctions=OUTPUT_DIR+"/J6/{sample}.J6.gz"
    params:
        entropy=config["entropy"],
        total_count=config["total_count"],
        gtag=("", "-g")[config["gtag"]]
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
    shell:
         "python3 -m workflow.scripts.merge_junctions "
         "{input.stranded_junctions} "
         "-o {output.merged_junctions}"