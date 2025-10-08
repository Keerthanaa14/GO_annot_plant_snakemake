rule annotate_diamond_uniprot_online:
    input:
        diamond="results/diamond/{sample}.tsv"
    output:
        annotated="results/diamond/{sample}.annotated.tsv"
    log:
        "logs/annotate_diamond_{sample}.log"
    run:
        import pandas as pd
        import requests, time, sys

        df = pd.read_csv(input.diamond, sep="\t", header=None)
        df.columns = ["qseqid","sseqid","pident","length","mismatch","gapopen",
                      "qstart","qend","sstart","send","evalue","bitscore","stitle"]

        # Query UniProt online (RefSeq -> UniProt)
        def get_uniprot_info(ncbi_ids):
            url = "https://rest.uniprot.org/idmapping/run"
            results = {}
            batch_size = 100
            for i in range(0, len(ncbi_ids), batch_size):
                batch = ncbi_ids[i:i+batch_size]
                data = {"from":"RefSeq_Protein","to":"UniProtKB","ids":",".join(batch)}
                r = requests.post(url, data=data)
                if r.status_code != 200:
                    print(f"Error {r.status_code}", file=sys.stderr)
                    continue
                job_id = r.json()["jobId"]
                status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
                while True:
                    status = requests.get(status_url).json()
                    if status.get("jobStatus") == "FINISHED":
                        break
                    time.sleep(1)
                results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}?format=tsv"
                res = requests.get(results_url)
                if res.status_code != 200: continue
                for line in res.text.strip().split("\n")[1:]:
                    from_id, to_id, db = line.split("\t")
                    results[from_id] = {"uniprot_id": to_id, "db_source": db}
                time.sleep(1)
            return results

        ncbi_ids = df["sseqid"].unique().tolist()
        uniprot_map = get_uniprot_info(ncbi_ids)

        df["uniprot_id"] = df["sseqid"].map(lambda x: uniprot_map.get(x, {}).get("uniprot_id","NA"))
        df["db_source"] = df["sseqid"].map(lambda x: uniprot_map.get(x, {}).get("db_source","NA"))

        df.to_csv(output.annotated, sep="\t", index=False)
        with open(log[0], "a") as l:
            l.write(f"[ANNOTATE {wildcards.sample}] Completed\n")