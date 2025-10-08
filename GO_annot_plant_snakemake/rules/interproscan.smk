rule interproscan:
    input:
        fasta="results/fasta/{sample}_top.fasta",
        eggnog="results/eggnog/{sample}.tsv.emapper.annotations"
    output:
        tsv="results/interpro/{sample}.tsv"
    log:
        "logs/interpro_{sample}.log"
    params:
        opts=config.get("interproscan_opts", "")
    shell:
        """
        mkdir -p results/interpro results/interpro/tmp
        module load biokit
        module load interproscan
        echo "[INTERPROSCAN {wildcards.sample}] $(date)" >> {log}
        cluster_interproscan \
            -i {input.fasta} \
            -o {output.tsv} \
            -f TSV \
            -dp results/interpro/tmp/{wildcards.sample} \
            {params.opts} \
            &>> {log}
        """