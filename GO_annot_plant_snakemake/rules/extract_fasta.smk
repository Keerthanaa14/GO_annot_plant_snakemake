rule extract_fasta:
    input:
        diamond_out="results/diamond/{sample}.tsv",
        fasta="/projappl/project_2015371/nrplants.fasta"
    output:
        "results/fasta/{sample}_top.fasta"
    params:
        ids="results/fasta/{sample}_top.ids"
    shell:
        """
        mkdir -p results/fasta
        cut -f2 {input.diamond_out} | sort -u > {params.ids}
        seqkit grep -f {params.ids} {input.fasta} > {output}
        """