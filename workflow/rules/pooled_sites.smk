rule count_pooled_sites:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=rules.index_bam.output.bam_index,
        junctions=rules.merge_junctions.output.merged_junctions,
        stats=rules.count_junctions.output.library_stats
    output:
        pooled_sites=OUTPUT_DIR+"/PS1/{sample}.PS1.gz"
    params:
        threads=THREADS,
        primary=("", "-p")[config["primary"]],
        unique=("", "-u")[config["unique"]]
    shell:
        "python3 -m workflow.scripts.count_sites "
        "-i {input.bam} "
        "-j {input.junctions} "
        "-s {input.stats} "
        "-o {output.pooled_sites} "
        "{params.primary} {params.unique} "
        "-t {params.threads}"


rule aggregate_pooled_sites:
    input:
        sites=rules.count_pooled_sites.output.pooled_sites,
        stats=rules.count_pooled_sites.input.stats
    output:
        aggregated_pooled_sites=OUTPUT_DIR+"/PS2/{sample}.PS2.gz"
    params:
        min_offset=config["min_offset"]
    shell:
        "python3 -m workflow.scripts.aggregate_sites "
        "-i {input.sites} "
        "-s {input.stats} "
        "-o {output.aggregated_pooled_sites} "
        "-m {params.min_offset}"


rule filter_pooled_sites:
    input:
        aggregated_pooled_sites=rules.aggregate_pooled_sites.output.aggregated_pooled_sites
    output:
        filtered_pooled_sites=OUTPUT_DIR+"/PS6/{sample}.PS6.gz"
    params:
        entropy=config["entropy"],
        total_count=config["total_count"]
    shell:
         "python3 -m workflow.scripts.filter "
         "-i {input.aggregated_pooled_sites} "
         "--sites "
         "-e {params.entropy} "
         "-c {params.total_count} "
         "-o {output.filtered_pooled_sites}"
