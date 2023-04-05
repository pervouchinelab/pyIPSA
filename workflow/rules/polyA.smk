rule count_polyA:
    input:
        bam=INPUT_DIR+"/{sample}.bam",
        bam_index=rules.index_bam.output.bam_index,
        library_stats=rules.count_junctions.output.library_stats
    output:
        polyA=OUTPUT_DIR+"/A1/{sample}.A1.gz",
    threads: THREADS
    params:
        primary=("", "-p")[config["primary"]],
        unique=("", "-u")[config["unique"]]
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.count_polyA "
        "-i {input.bam} "
        "-o {output.polyA} "
        "{params.primary} {params.unique} "
        "-t {threads}"


rule aggregate_polyA:
    input:
        polyA=rules.count_polyA.output.polyA,
        library_stats=rules.count_junctions.output.library_stats
    output:
        aggregated_polyA=OUTPUT_DIR+"/A2/{sample}.A2.gz"
    params:
        min_overhang=config["min_overhang"],
    conda: "../envs/scripts-common.yaml"
    shell:
        "python3 -m workflow.scripts.aggregate_polyA "
        "-i {input.polyA} "
        "-s {input.library_stats} "
        "-o {output.aggregated_polyA} "
        "--min_overhang {params.min_overhang} "