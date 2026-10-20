# host_comm_meta.smk
# Split metatranscriptome reads (fastp output of meta_transcriptome.smk) into
# host reads and community (non-host) reads, same logic as host_community.smk.
#
# Standalone run:  snakemake -s host_comm_meta.smk -n   (dry run)
#                  snakemake -s host_comm_meta.smk --cores 8
#
# Input : {MT_OUT_DIR}/fastp/{SRA_ID}_{1,2}.fastq.gz  (from meta_transcriptome.smk)
#         {ref_dir}/{taxid}.fna                       (from download_genome.smk)
# Output: {MT_OUT_DIR}/host_community/{SRA_ID}_host_{1,2}.fastq.gz
#         {MT_OUT_DIR}/host_community/{SRA_ID}_community_{1,2}.fastq.gz
#
# NOTE: `rule all` below is for standalone runs. If this file is included from
# the root Snakefile, delete it and add MT_TARGETS to the root `rule all`.

configfile: "config/config.yaml"

# MT_* prefix: .smk files share one global namespace, avoid clobbering
# OUT_DIR / INP_DIR / SRA_ID defined in the other rule files.
MT_SRA_IDS = config["metatranscriptome"]["sra"]["sra_id"]
MT_OUT_DIR = config["metatranscriptome"]["output_dir"]
MT_FASTP_DIR = f"{MT_OUT_DIR}/fastp"
MT_SPLIT_DIR = f"{MT_OUT_DIR}/host_community"

MT_TAXID = config["metatranscriptome"]["taxid"]
MT_REFERENCE = config["metatranscriptome"]["ref_dir"]
MT_HOST_GENOME = f"{MT_REFERENCE}/{MT_TAXID}.fna"
MT_BWA_INDEX = expand(
    "{reference}/{taxid}.fna.{ext}",
    reference=MT_REFERENCE,
    taxid=MT_TAXID,
    ext=["amb", "ann", "bwt", "pac", "sa"],
)

MT_TARGETS = expand(
    MT_SPLIT_DIR + "/{SRA_ID}_{part}_{read}.fastq.gz",
    SRA_ID=MT_SRA_IDS,
    part=["host", "community"],
    read=[1, 2],
)


rule all:
    input:
        MT_TARGETS


# bwa index of the host genome (writes <genome>.amb/.ann/.bwt/.pac/.sa)
rule bwa_index_metatranscriptome:
    input:
        host_genome = MT_HOST_GENOME
    output:
        MT_BWA_INDEX
    log:
        f"{MT_OUT_DIR}/logs/bwa_index.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {MT_OUT_DIR}/logs
        bwa index {input.host_genome} 2> {log}
        """


# align (for each SRA ID) + sort + index
rule bwa_align_metatranscriptome:
    input:
        r1 = f"{MT_FASTP_DIR}/{{SRA_ID}}_1.fastq.gz",
        r2 = f"{MT_FASTP_DIR}/{{SRA_ID}}_2.fastq.gz",
        host_genome = MT_HOST_GENOME,
        index_files = MT_BWA_INDEX
    output:
        bam = f"{MT_SPLIT_DIR}/{{SRA_ID}}_aligned_sorted.bam",
        bai = f"{MT_SPLIT_DIR}/{{SRA_ID}}_aligned_sorted.bam.bai"
    threads:
        config["threads"]["bwa_align"]
    params:
        sort_threads = 4,
        sort_memory = "2G",
        extra_params = "-v 1"
    log:
        f"{MT_SPLIT_DIR}/logs/{{SRA_ID}}_bwa_align.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {MT_SPLIT_DIR}/logs

        # bwa mem -> sort without an intermediate SAM file
        bwa mem -t {threads} {params.extra_params} \\
            {input.host_genome} {input.r1} {input.r2} 2> {log} \\
        | samtools sort -@ {params.sort_threads} -m {params.sort_memory} \\
            -T {output.bam}.tmp -o {output.bam} - 2>> {log}

        samtools index -@ {params.sort_threads} {output.bam} 2>> {log}
        """


# split into host / community, keeping read pairs
rule split_reads_metatranscriptome:
    input:
        bam = f"{MT_SPLIT_DIR}/{{SRA_ID}}_aligned_sorted.bam",
        bai = f"{MT_SPLIT_DIR}/{{SRA_ID}}_aligned_sorted.bam.bai"
    output:
        host_r1 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_host_1.fastq.gz",
        host_r2 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_host_2.fastq.gz",
        community_r1 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_community_1.fastq.gz",
        community_r2 = f"{MT_SPLIT_DIR}/{{SRA_ID}}_community_2.fastq.gz"
    threads:
        config["threads"]["split_reads"]
    params:
        collated = f"{MT_SPLIT_DIR}/{{SRA_ID}}_collated.tmp.bam"
    log:
        f"{MT_SPLIT_DIR}/logs/{{SRA_ID}}_split_reads.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {MT_SPLIT_DIR}/logs

        # group mates together (samtools fastq needs name-collated input)
        samtools collate -@ {threads} -o {params.collated} \\
            {input.bam} {params.collated}.prefix 2> {log}

        # community: both mates unmapped (flag 4 + 8 = 12)
        samtools fastq -@ {threads} -n -f 12 \\
            -1 {output.community_r1} -2 {output.community_r2} \\
            -0 /dev/null -s /dev/null {params.collated} 2>> {log}

        # host: everything except pairs where both mates are unmapped
        samtools fastq -@ {threads} -n -G 12 \\
            -1 {output.host_r1} -2 {output.host_r2} \\
            -0 /dev/null -s /dev/null {params.collated} 2>> {log}

        rm -f {params.collated}
        """

