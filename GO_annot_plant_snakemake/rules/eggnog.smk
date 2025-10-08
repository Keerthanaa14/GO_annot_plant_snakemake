rule eggnog_mapper:
    input:
        fasta="results/fasta/{sample}_top.fasta"
    output:
        "results/eggnog/{sample}.tsv.emapper.annotations"
    log:
        "logs/eggnog_{sample}.log"
    params:
        emapper=config["eggnog_mapper_path"],
        data_dir=config["eggnog_db"],
        opts=config["eggnog_opts"]
    shell:
        """
        mkdir -p results/eggnog
        echo "[EGGNOG {wildcards.sample}] $(date)" >> {log}
        source /projappl/project_2015371/eggnog-mapper/venv/bin/activate
        python {params.emapper} \
            -i {input.fasta} \
            --data_dir {params.data_dir} \
            {params.opts} \
            -o results/eggnog/{wildcards.sample}.tsv \
            &>> {log}
        """