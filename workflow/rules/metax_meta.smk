# metax_metatranscriptome.smk
# Taxonomic profiling of the community (non-host) metatranscriptome reads with
# Metax (https://github.com/hzi-bifo/Metax), mode: precision.
#
# Standalone run:  snakemake -s metax_metatranscriptome.smk -n   (dry run)
#                  snakemake -s metax_metatranscriptome.smk --cores 8
#
# Input : {MT_OUT_DIR}/host_community/{SRA_ID}_community_{1,2}.fastq.gz
#         (from host_comm_meta.smk)
# Output: {MT_OUT_DIR}/metax/{SRA_ID}/{SRA_ID}.profile.txt
#         (Metax also writes .sam.gz, .classify.txt and .log next to it)
#
# Metax runs in its own conda env (default name: `metax`) via `conda run`, so
# this file can be launched from the `snakemake` env. Create it once:
#   conda create -n metax -c zldeng metax
#
# config.yaml (under `metatranscriptome:`):
#   metax:
#     db_json: "/path/to/metax_db.json"
#     dmp_dir: "/path/to/taxdump"      # nodes.dmp, names.dmp, merged.dmp
#     mode: "precision"
#     sequencer: "Illumina"
#     threads: 8
#     conda_env: "metax"               # optional, default "metax"
#
# NOTE: `rule all` below is for standalone runs. If this file is included from
# the root Snakefile, delete it and add MT_METAX_TARGETS to the root `rule all`.

configfile: "config/config.yaml"

MT_SRA_IDS = config["metatranscriptome"]["sra"]["sra_id"]
MT_OUT_DIR = config["metatranscriptome"]["output_dir"]
MT_SPLIT_DIR = f"{MT_OUT_DIR}/host_community"
MT_METAX_DIR = f"{MT_OUT_DIR}/metax"

MT_METAX = config["metatranscriptome"]["metax"]
MT_METAX_DB = MT_METAX["db_json"]
MT_METAX_DMP_DIR = MT_METAX["dmp_dir"]
MT_METAX_ENV = MT_METAX.get("conda_env", "metax")
MT_METAX_DMP_FILES = expand(
    "{dmp_dir}/{name}.dmp",
    dmp_dir=MT_METAX_DMP_DIR,
    name=["nodes", "names", "merged"],
)

MT_METAX_TARGETS = expand(
    MT_METAX_DIR + "/{SRA_ID}/{SRA_ID}.profile.txt",
    SRA_ID=MT_SRA_IDS,
)


rule all:
    input:
        MT_METAX_TARGETS


rule metax_profile_metatranscriptome:
    input:
        r1 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_community_1.fastq.gz",
        r2 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_community_2.fastq.gz",
        db = MT_METAX_DB,
        dmp = MT_METAX_DMP_FILES
    output:
        profile = f"{MT_METAX_DIR}/{{SRA_ID}}/{{SRA_ID}}.profile.txt"
    threads:
        MT_METAX["threads"]
    params:
        outdir = f"{MT_METAX_DIR}/{{SRA_ID}}",
        prefix = f"{MT_METAX_DIR}/{{SRA_ID}}/{{SRA_ID}}",
        dmp_dir = MT_METAX_DMP_DIR,
        mode = MT_METAX["mode"],
        sequencer = MT_METAX["sequencer"],
        conda_env = MT_METAX_ENV
    log:
        f"{MT_METAX_DIR}/logs/{{SRA_ID}}_metax.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir} {MT_METAX_DIR}/logs

        conda run -n {params.conda_env} --no-capture-output metax profile \\
            --db {input.db} \\
            --dmp-dir {params.dmp_dir} \\
            --in-seq {input.r1},{input.r2} \\
            --is-paired \\
            --sequencer {params.sequencer} \\
            --mode {params.mode} \\
            --outprefix {params.prefix} \\
            --threads {threads} \\
            --compress-sam > {log} 2>&1
        """
