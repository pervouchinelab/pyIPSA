rule count_sites:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=rules.index_bam.output.bam_index,
        junctions=rules.choose_strand.output.stranded_junctions,
        stats=rules.count_junctions.output.library_stats
    output:
        sites=OUTPUT_DIR+"/S1/{sample}.S1.gz"
    threads: THREADS
    params:
        primary=("", "-p")[config["primary"]],
        unique=("", "-u")[config["unique"]]
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.count_sites "
        "-i {input.bam} "
        "-j {input.junctions} "
        "-s {input.stats} "
        "-o {output.sites} "
        "{params.primary} {params.unique} "
        "-t {threads}"


rule aggregate_sites:
    input:
        sites=rules.count_sites.output.sites,
        stats=rules.count_sites.input.stats
    output:
        aggregated_sites=OUTPUT_DIR+"/S2/{sample}.S2.gz"
    params:
        min_offset=config["min_offset"]
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.aggregate_sites "
        "-i {input.sites} "
        "-s {input.stats} "
        "-o {output.aggregated_sites} "
        "-m {params.min_offset}"


rule filter_sites:
    input:
        aggregated_sites=rules.aggregate_sites.output.aggregated_sites
    output:
        filtered_sites=OUTPUT_DIR+"/S6/{sample}.S6.gz"
    params:
        entropy=config["entropy"],
        total_count=config["total_count"]
    conda: "../envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.filter "
         "-i {input.aggregated_sites} "
         "--sites "
         "-e {params.entropy} "
         "-c {params.total_count} "
         "-o {output.filtered_sites}"