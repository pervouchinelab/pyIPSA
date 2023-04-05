rule compute_rates:
    input:
        filtered_junctions=rules.filter_junctions.output.filtered_junctions,
        filtered_sites=rules.filter_sites.output.filtered_sites
    output:
        rates=OUTPUT_DIR+"/R/{sample}.R.gz"
    conda: "envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.compute_rates "
         "-j {input.filtered_junctions} "
         "-s {input.filtered_sites} "
         "-o {output.rates}"


rule compute_pooled_rates:
    input:
        filtered_junctions=rules.filter_junctions.output.filtered_junctions,
        filtered_pooled_sites=rules.filter_pooled_sites.output.filtered_pooled_sites
    output:
        rates=OUTPUT_DIR+"/PR/{sample}.PR.gz"
    conda: "envs/scripts-common.yaml"
    shell:
         "python3 -m workflow.scripts.compute_rates "
         "-j {input.filtered_junctions} "
         "-s {input.filtered_pooled_sites} "
         "-o {output.rates}"