import os, glob, pathlib, shutil, csv
from collections import defaultdict

configfile: "config.yaml"

OUT       = config["out_dir"]
MAGS_DIR  = config["mags_dir"]
THREADS   = int(config.get("threads", 8))
DBS       = config["dbs"]
PROKKA    = config.get("prokka", {})
OVERRIDES = config.get("kingdom_overrides", {})

# Find MAGs (fa/fna/fasta)
MAG_FASTAS = sorted(
    glob.glob(os.path.join(MAGS_DIR, "*.fa")) +
    glob.glob(os.path.join(MAGS_DIR, "*.fna")) +
    glob.glob(os.path.join(MAGS_DIR, "*.fasta"))
)
assert MAG_FASTAS, f"No MAG FASTA files found under {MAGS_DIR}"

MAG_IDS = [pathlib.Path(p).stem for p in MAG_FASTAS]
MAG_MAP = {pathlib.Path(p).stem: p for p in MAG_FASTAS}

def prokka_outputs(m):
    base = f"{OUT}/{m}/prokka"
    return [
        f"{base}/{m}.faa",   # proteins
        f"{base}/{m}.gff",   # features
        f"{base}/{m}.ffn",   # CDS nucleotides
    ]

def emapper_outputs(m):
    base = f"{OUT}/{m}/eggnog"
    return [
        f"{base}/{m}.emapper.annotations.tsv",
        f"{base}/{m}.emapper.seed_orthologs",
    ]

def summary_outputs(m):
    base = f"{OUT}/{m}/summary"
    return [
        f"{base}/{m}.KO_counts.tsv",
        f"{base}/{m}.COG_category_counts.tsv",
    ]

rule all:
    input:
        # per-MAG annotations and summaries
        sum([prokka_outputs(m) for m in MAG_IDS], []) +
        sum([emapper_outputs(m) for m in MAG_IDS], []) +
        sum([summary_outputs(m) for m in MAG_IDS], []) + [
            # combined cohort summaries
            f"{OUT}/combined/KO_counts_wide.tsv",
            f"{OUT}/combined/COG_category_counts_wide.tsv",
        ]

# -------------------------
# Gene prediction (Prokka)
# -------------------------
rule prokka:
    input:
        fasta=lambda wc: MAG_MAP[wc.mag]
    output:
        faa = f"{OUT}/{{mag}}/prokka/{{mag}}.faa",
        gff = f"{OUT}/{{mag}}/prokka/{{mag}}.gff",
        ffn = f"{OUT}/{{mag}}/prokka/{{mag}}.ffn"
    conda: "envs/prokka.yaml"
    threads: 8
    resources:
        mem_mb=16000
    run:
        os.makedirs(f"{OUT}/{wildcards.mag}/prokka", exist_ok=True)
        kingdom = OVERRIDES.get(wildcards.mag, PROKKA.get("kingdom_default", "Bacteria"))
        gram    = PROKKA.get("gram", "")
        genus   = PROKKA.get("genus", "")
        species = PROKKA.get("species", "")
        gram_opt  = f"--gram {gram}" if gram else ""
        genus_opt = f"--genus {genus}" if genus else ""
        species_opt = f"--species {species}" if species else ""

        shell("""
            prokka --outdir {OUT}/{wildcards.mag}/prokka --prefix {wildcards.mag} \
                   --kingdom {kingdom} --cpus {threads} \
                   {gram_opt} {genus_opt} {species_opt} \
                   {input.fasta}
        """)

# -------------------------
# Functional annotation (eggNOG-mapper)
# -------------------------
rule eggnog_mapper:
    input:
        faa = f"{OUT}/{{mag}}/prokka/{{mag}}.faa"
    output:
        tsv   = f"{OUT}/{{mag}}/eggnog/{{mag}}.emapper.annotations.tsv",
        seeds = f"{OUT}/{{mag}}/eggnog/{{mag}}.emapper.seed_orthologs"
    conda: "envs/emapper.yaml"
    threads: 8
    shell:
        r"""
        mkdir -p {OUT}/{wildcards.mag}/eggnog
        if [[ -s {input.faa} ]]; then
            emapper.py -i {input.faa} -o {wildcards.mag}.emapper \
                --output_dir {OUT}/{wildcards.mag}/eggnog \
                --data_dir {DBS[eggnog_data_dir]} \
                --cpu {threads} --override \
                --decorate_gff no \
                --target_orthologs {config[eggnog][target_orthologs]} \
                --tax_scope {config[eggnog][tax_scope]} \
                --go_evidence {config[eggnog][go_evidence]} \
                --usemem \
                {"--m " + config["eggnog"]["mode"] if config["eggnog"].get("mode") else ""}
        else
            # empty inputs -> create empty outputs
            : > {output.tsv}
            : > {output.seeds}
        fi
        """

# -------------------------
# Per-MAG summaries
# -------------------------
rule summarize_emapper_per_mag:
    input:
        tsv = f"{OUT}/{{mag}}/eggnog/{{mag}}.emapper.annotations.tsv"
    output:
        ko = f"{OUT}/{{mag}}/summary/{{mag}}.KO_counts.tsv",
        cog = f"{OUT}/{{mag}}/summary/{{mag}}.COG_category_counts.tsv"
    conda: "envs/emapper.yaml"
    run:
        import pandas as pd
        os.makedirs(f"{OUT}/{wildcards.mag}/summary", exist_ok=True)
        if not os.path.getsize(input.tsv):
            # empty placeholders
            pd.DataFrame(columns=["KO","count"]).to_csv(output.ko, sep="\t", index=False)
            pd.DataFrame(columns=["COG_category","count"]).to_csv(output.cog, sep="\t", index=False)
        else:
            # eggNOG tsv has columns incl. 'KEGG_ko' and 'COG_category' (comma/pipe-separated)
            df = pd.read_csv(input.tsv, sep="\t", comment="#", low_memory=False)
            # KEGG KO counts
            if "KEGG_ko" in df.columns:
                kos = (
                    df["KEGG_ko"]
                    .dropna()
                    .astype(str)
                    .str.replace(" ", "", regex=False)
                    .str.split(",|\\|", regex=True)
                    .explode()
                )
                ko_counts = kos[kos.str.len() > 0].value_counts().rename_axis("KO").reset_index(name="count")
            else:
                ko_counts = pd.DataFrame(columns=["KO","count"])
            ko_counts.to_csv(output.ko, sep="\t", index=False)

            # COG category counts (single-letter categories, possibly multiple per gene)
            if "COG_category" in df.columns:
                cogs = (
                    df["COG_category"]
                    .dropna()
                    .astype(str)
                    .str.replace(" ", "", regex=False)
                    .str.split(",|\\|", regex=True)
                    .explode()
                )
                cog_counts = cogs[cogs.str.len() > 0].value_counts().rename_axis("COG_category").reset_index(name="count")
            else:
                cog_counts = pd.DataFrame(columns=["COG_category","count"])
            cog_counts.to_csv(output.cog, sep="\t", index=False)

# -------------------------
# Combined summaries (wide)
# -------------------------
rule combine_ko_counts:
    input:
        expand(f"{OUT}/{{mag}}/summary/{{mag}}.KO_counts.tsv", mag=MAG_IDS)
    output:
        wide = f"{OUT}/combined/KO_counts_wide.tsv"
    conda: "envs/emapper.yaml"
    run:
        import pandas as pd
        os.makedirs(f"{OUT}/combined", exist_ok=True)
        tables = []
        for mag in MAG_IDS:
            p = f"{OUT}/{mag}/summary/{mag}.KO_counts.tsv"
            df = pd.read_csv(p, sep="\t")
            if len(df):
                df = df.set_index("KO").rename(columns={"count": mag})
            else:
                df = pd.DataFrame(columns=[mag], index=[])
            tables.append(df)
        if tables:
            wide = pd.concat(tables, axis=1).fillna(0).astype(int)
        else:
            wide = pd.DataFrame()
        wide.to_csv(output.wide, sep="\t")

rule combine_cog_counts:
    input:
        expand(f"{OUT}/{{mag}}/summary/{{mag}}.COG_category_counts.tsv", mag=MAG_IDS)
    output:
        wide = f"{OUT}/combined/COG_category_counts_wide.tsv"
    conda: "envs/emapper.yaml"
    run:
        import pandas as pd
        os.makedirs(f"{OUT}/combined", exist_ok=True)
        tables = []
        for mag in MAG_IDS:
            p = f"{OUT}/{mag}/summary/{mag}.COG_category_counts.tsv"
            df = pd.read_csv(p, sep="\t")
            if len(df):
                df = df.set_index("COG_category").rename(columns={"count": mag})
            else:
                df = pd.DataFrame(columns=[mag], index=[])
            tables.append(df)
        if tables:
            wide = pd.concat(tables, axis=1).fillna(0).astype(int)
        else:
            wide = pd.DataFrame()
        wide.to_csv(output.wide, sep="\t")

