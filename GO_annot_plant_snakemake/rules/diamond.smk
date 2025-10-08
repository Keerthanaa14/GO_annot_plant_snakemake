rule diamond_blastp:
    input:
        fasta=lambda wc: config["samples"][wc.sample]
    output:
        "results/diamond/{sample}.tsv"
    log:
        "logs/diamond_{sample}.log"
    threads: config["diamond_threads"]
    params:
        db=config["diamond_db"],
        evalue=config["diamond_evalue"],
        top_hits=config["top_hits"],
        block_size=config["diamond_block"],
        chunks=config["diamond_chunks"]
    shell:
        """
        module load biokit
        module load diamond
        echo "[DIAMOND {wildcards.sample}] $(date)" >> {log}

        diamond blastp \
            --query {input.fasta} \
            --db {params.db} \
            --evalue {params.evalue} \
            --max-target-seqs {params.top_hits} \
            --threads {threads} \
            --block-size {params.block_size} \
            --index-chunks {params.chunks} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --out {output} \
            &>> {log}
        """