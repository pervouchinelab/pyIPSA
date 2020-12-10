rule download_genome:
    output:
        genome="genomes/{org}.fa"
    params:
        url=lambda wildcards: config["genome_urls"][wildcards.org]
    shell:
        """
        wget -O {output.genome}.gz {params.url}
        gunzip {output.genome}.gz
        """