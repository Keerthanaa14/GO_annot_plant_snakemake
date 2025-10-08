rule combine_results:
    input:
        diamond="results/diamond/{sample}.tsv",
        eggnog="results/eggnog/{sample}.tsv.emapper.annotations",
        interpro="results/interpro/{sample}.tsv",
        go_obo="data/go-basic.obo"
    output:
        "results/combined/{sample}_summary.tsv"
    log:
        "logs/combine_{sample}.log"
    run:
        import pandas as pd
        from pathlib import Path

        # ----------------------------
        # Create output folder
        # ----------------------------
        Path("results/combined").mkdir(parents=True, exist_ok=True)

        # ----------------------------
        # Helper: parse GO OBO → dict
        # ----------------------------
        def parse_go_obo(obo_file):
            go_map = {}
            go_id, go_name = None, None
            with open(obo_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("id: GO:"):
                        go_id = line.split("id: ")[1]
                    elif line.startswith("name: ") and go_id:
                        go_name = line.split("name: ")[1]
                    elif line == "" and go_id and go_name:
                        go_map[go_id] = go_name
                        go_id, go_name = None, None
                if go_id and go_name:
                    go_map[go_id] = go_name
            return go_map

        with open(log[0], "a") as l:
            l.write(f"[COMBINE {wildcards.sample}] Started at {pd.Timestamp.now()}\n")

            # ----------------------------
            # Load DIAMOND
            # ----------------------------
            diamond_cols = [
                "qseqid","sseqid","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore","stitle"
            ]
            d = pd.read_csv(input.diamond, sep="\t", header=None, names=diamond_cols)
            d["sseqid"] = d["sseqid"].astype(str).str.strip()
            d["qseqid"] = d["qseqid"].astype(str).str.strip()
            d["sseqid_simple"] = d["sseqid"].str.split("|").str[-1]
            l.write(f"DIAMOND loaded: {len(d)} rows\n")

            # ----------------------------
            # Load EggNOG
            # ----------------------------
            eggnog_cols = [
                "query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl",
                "COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",
                "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE",
                "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"
            ]
            e = pd.read_csv(
                input.eggnog,
                sep="\t",
                comment="#",
                header=None,
                names=eggnog_cols,
                low_memory=False
            )
            e["query"] = e["query"].astype(str).str.strip()
            e["query_simple"] = e["query"].str.split("|").str[-1]
            l.write(f"EggNOG loaded: {len(e)} rows\n")

            # ----------------------------
            # Load InterProScan
            # ----------------------------
            interpro_cols = [
                "protein_id","md5","length","analysis","signature_acc","signature_desc",
                "start","end","score","status","date","ipr_acc","ipr_desc","go_terms","pathways"
            ]
            i = pd.read_csv(input.interpro, sep="\t", header=None, names=interpro_cols, dtype=str)
            i["protein_id"] = i["protein_id"].astype(str).str.strip()
            i["protein_id_simple"] = i["protein_id"].str.split("|").str[-1]
            l.write(f"InterProScan loaded: {len(i)} rows\n")

            # ----------------------------
            # Collapse InterProScan per protein
            # ----------------------------
            i_collapsed = i.groupby("protein_id_simple").agg({
                "ipr_desc": lambda x: ";".join(sorted(set(x.dropna()))),
                "go_terms": lambda x: ";".join(sorted(set(x.dropna()))),
                "pathways": lambda x: ";".join(sorted(set(x.dropna())))
            }).reset_index()
            l.write(f"InterProScan collapsed: {len(i_collapsed)} proteins\n")

            # ----------------------------
            # Parse GO mapping from OBO
            # ----------------------------
            go_map = parse_go_obo(input.go_obo)
            l.write(f"GO mapping parsed: {len(go_map)} terms\n")

            # ----------------------------
            # Map GO IDs → descriptions
            # ----------------------------
            def map_go_desc_interpro(go_str):
                if pd.isna(go_str) or go_str == "":
                    return ""
                terms = go_str.replace("|", ";").split(";")
                descs = []
                for t in terms:
                    t_clean = t.split("(")[0].strip()
                    descs.append(go_map.get(t_clean, t_clean))
                return ";".join(sorted(set(descs)))

            def map_go_desc_eggnog(go_str):
                if pd.isna(go_str) or go_str == "":
                    return ""
                terms = go_str.replace(",", ";").replace("|", ";").split(";")
                descs = []
                for t in terms:
                    t_clean = t.strip()
                    descs.append(go_map.get(t_clean, t_clean))
                return ";".join(sorted(set(descs)))

            # Optionally, create mapped GO columns
            i_collapsed["go_desc"] = i_collapsed["go_terms"].apply(map_go_desc_interpro)
            e["go_desc"] = e["GOs"].apply(map_go_desc_eggnog)

            # ----------------------------
            # Merge DIAMOND → InterProScan → EggNOG
            # ----------------------------
            merged = d.merge(
                i_collapsed, left_on="sseqid_simple", right_on="protein_id_simple", how="left"
            ).merge(
                e, left_on="sseqid_simple", right_on="query_simple", how="left"
            )
            l.write(f"Merged table: {len(merged)} rows\n")

            # ----------------------------
            # Save output
            # ----------------------------
            merged.to_csv(output[0], sep="\t", index=False)
            l.write(f"[COMBINE {wildcards.sample}] Finished at {pd.Timestamp.now()}\n")