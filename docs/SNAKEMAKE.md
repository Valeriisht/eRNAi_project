# Snakemake part review

Here we look at the `.smk` files in `rules`.

## A bit about Snakemake

| Construct | Meaning |
| --- | --- |
| `configfile:` | read a YAML file and make its values available as `config[...]` |
| `rule NAME:` | name a computation step |
| `input:` | files that must exist before the step runs |
| `output:` | files whose appearance means the step completed successfully |
| `wildcard` | a placeholder, e.g. `{sra_id}` becomes `SRR8265535` |
| `shell:` | the command Snakemake will pass to the shell |

Snakemake builds its graph not by line order but by matching input and output names. If `rule B` requests a file produced by `rule A`, B waits for A.

## `config/config.yaml`

### Paths

- `output_dir: "kraken"` — the directory the rules expect to write results into. Not only Kraken files end up there.
- `input_dir: "metagenome/dehost"` — the directory of expected reads after host sequences have been removed.
- `ref_dir` — the directory holding the host reference.
- `taxid: "galGal6"` — the assembly/reference name.
- `sample_name` — the shared prefix for the merged report.

### SRA and resources

- `sra.sra_id` — the list of metagenomic SRA runs.
- `paired: True` means paired-end reads are expected.
- `thread: 2` — threads for `fasterq-dump`.

### BWA/SAMtools resources

The `threads:` block is declared twice. YAML will keep the **second** block, so the values actually in effect are 8, 8, 8, 4. This should be collapsed into a single block so the true settings are documented.

### Kallisto

- `kmer_size: 31` — the desired k-mer length for the Kallisto index.
- `bootstrap: 30` — the number of bootstrap replicates used to estimate uncertainty in the abundance estimate.

### Metagenome/QC

- `algorithm` selects the Kraken2 or MetaPhlAn branch (NOT IMPLEMENTED).
- `database` — an absolute path — SHOULD BE REMOVED.
- `read_length: 150` — the read length passed to Bracken.
- `taxonomic_level: ['S']` — species level is requested. The code also supports `G` and `P`, but the current configuration doesn't produce them.
- fastp thresholds: Q20 and 50 nt. `detect_adapters` is set, but the rule always enables adapter detection and never checks this boolean.

---

## `Snakefile`: what happens here

### Main script

`configfile: "config/config.yaml"` makes the YAML values available as the `config` dictionary.

### Global variables

Each line copies a convenient value out of `config`:

```python
INP_DIR  # path to input reads
OUT_DIR  # path to results
SAMPLE   # name of the merged report
SRA_IDS  # list of accession IDs
DB       # Kraken database
ALGO     # kraken2/metaphlan
READ_LEN # read length for Bracken
LEVELS   # S/G/P
```

### Metagenome target

- Includes the rules from `meta_genome.smk`.
- Declares `rule all`, i.e. the desired final files.
- `expand(...)` multiplies the template across all SRA IDs — for two IDs, for example, it produces two Kraken report names.

### Host removal branch

- Includes `host_community.smk`.
- The template produces two files per SRA ID: host reads and metagenome reads.

### QC and transcriptome

- Includes two modules.

### Download genome

- The rule for downloading genomes.

### Honest assessment of the `Snakefile`

This is a draft file that bundles several alternative runs into a single document.

The correct structure is to run things in this order (they were previously run separately):

```text
Snakefile
  ├── include: preprocessing.smk
  ├── include: host_removal.smk
  ├── include: metagenome.smk
  ├── include: host_transcriptome.smk
  ├── include: metatranscriptome.smk
  └── one rule all: all selected final outputs
```

---

Now a walkthrough of the individual files — a bit about each one:

| Rule / block | Why it's needed |
| --- | --- |
| `rule all` | Defines the final outputs of the whole workflow and triggers the build of all dependencies. |
| `download_genome.smk` | Downloads the host reference genome/transcriptome for alignment and quantification. |
| `download_metagenome` | Downloads the raw metagenome FASTQ files by accession. |
| `pre_prossecing.smk` | Cleans reads: adapters, low quality, short sequences; produces QC reports. |
| Host genome indexing | Prepares the reference for fast alignment of metagenomic reads to the host. |
| `host_community.smk`: host-read removal | Removes host reads from the metagenome, leaving the potentially microbial reads. |
| `meta_genome.smk`: Kraken2 | Determines the taxonomic composition of microorganisms in the sample. |
| `meta_genome.smk`: Bracken | Estimates taxon abundance from the Kraken2 results. |
| Merging metagenome profiles | Builds a "sample × microorganism" matrix for integrative analysis. |
| `transcriptome.smk`: Kallisto index | Builds the host transcriptome index. |
| `transcriptome.smk`: Kallisto quantification | Estimates host transcript/gene expression from RNA-seq. |
| Merging RNA-seq results | Builds a "sample × host gene" matrix. |

## `workflow/rules/download_genome.smk`

### Parameters

- The config is re-read here. Under `include`, this is redundant but usually not critical.
- The reference ID and output directory are copied.
- The final filename `kraken/galGal6.fna` is computed.

### `download_genome`

- Declares the rule.
- Declares the named output `zip`.
- Passes `TAXID` into the command as `params.taxid`.
- Points to the log file.
- Runs the NCBI Datasets CLI and redirects stdout/stderr into the log. `--reference` requests the reference genome.

### `extract_genome`

- References the named output of the previous rule; this is a graph edge.
- `directory(...)` tells Snakemake the output is a directory, not a single file.
- Unzips the archive into that directory.

### `find_rename_genome`

- Expects the directory from `extract_genome`.
- Sets the final `.fna` target.
- Looks for the **first** `.fna` in the archive; if there are several, the choice is undocumented and could be wrong.
- Fails the job with an error if no fasta is found.
- Moves the found fasta to the standard output name.

### `clean_temp`

- Once `GENOME_FILE` exists, a `clean.done` marker is created.
- The command removes the zip and the unpacked directory.
- This is only safe if `GENOME_FILE` genuinely already holds the needed copy of the fasta. The documentation should note that intermediate files are deleted.

---

## `workflow/rules/pre_prossecing.smk`

### Setup

- Reads the shared config.
- Extracts the output directory and the list of SRA IDs. The name `SRA_ID` is misleading: it holds a list, not a single ID.
- Creates the output and `logs` directories **while the Snakefile is being parsed**, not inside a job. This is a side effect and should probably be removed in a future version — the output paths/commands themselves should create the directory.

### `prefetch_data`

- Writes `{SRA_ID}` in uppercase. Snakemake technically allows this as a wildcard name, but the rest of the project uses `{sra_id}`.
- Defines what gets substituted for the wildcard in the command.
- Downloads a single SRA archive.
- Has no input: this is a starting rule, triggered by the requested `.sra` output.

### `download_data`

- Connects this rule to `prefetch_data`: its input matches that rule's output.
- Names mate 1 and mate 2; the variable `r1` should actually be named `r2` — as written it makes the code confusing to read.
- Takes the SRA ID, thread count, and expected layout.
- Enables strict shell mode: an error, an uninitialized variable, or a pipe failure stops the command.
- Converts `.sra` to FASTQ. `strace` isn't part of the actual processing and requires Linux; it was likely left in for debugging and should be removed for portability.
- Renames the FASTQ files to the expected Snakemake output names.

### `process_paired_data`

- Takes the raw pair from the previous rule.
- Defines the two filtered FASTQ files and the JSON QC report.
- Passes the fastp parameters from the YAML. `detect_adapters` is declared here but isn't actually used as a condition in the shell command.
- Runs fastp: `-i/-I` are the two inputs, `-o/-O` the two outputs, `-q` the quality cutoff, `--length_required` the minimum length.
- Deletes the raw reads. This saves space but makes reprocessing and auditing harder; using `temp()` or storing the raw files separately would be better.

---

## `workflow/rules/host_community.smk`

- Dehosting.

### Parameters and side effect

- Copies the config. `SRA_ID` again holds a list and isn't used in this file.
- Creates the logs folder at parse time; see the note above.

### `bwa_index`

- Expects the host fasta `ref_dir/galGal6.fa`.
- Uses `expand()` for the five BWA index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) — a correct way to represent a single multi-part index.
- Pulls a resource value from the YAML, but the command **doesn't actually pass** the thread count: `bwa index` in this form doesn't take that parameter.
- Builds the index.

### `bwa_align`

- Takes the QC-cleaned paired reads. Their directory should match the preprocessing `output_dir`, or there should be an explicit move step; currently the YAML specifies different paths.
- Makes the index dependency explicit.
- Produces a SAM file for each `{sra_id}`.
- Sets `-v 1`, but that's a parameter of the old `bwa aln`, not `bwa mem`; it should be removed or replaced with a meaningful `bwa mem` option.
- Runs paired-end alignment: `-t` for threads, then the reference, R1, R2; the SAM goes to output, diagnostic stderr to the log.

### `sam_to_bam`

- Input is a SAM file; output is a coordinate-sorted BAM and its index.
- Hard-coded resources; better moved into the config.
- `samtools view` converts SAM into a temporary BAM.
- Sorts the BAM.
- Creates the `.bai` index.
- Removes the temporary BAM.
- `|| exit 1` after each operation prevents continuing on failure — this is good practice.

### `split_reads`

Idea: split reads by their alignment flag.

- `-f 4` selects reads with the `unmapped` flag set.
- `-F 4` excludes unmapped reads, i.e. keeps mapped ones.
- Writes unmapped mate 1 out as `metagenome_reads`.
- Writes mapped mate 1 out as `host_reads`.
- In both calls, `-2 /dev/null` **discards mate 2**. This can turn paired-end data into single-end data and break downstream rules.
- The outputs are named `.fastq.gz`, but `samtools fastq` isn't given any compression option and doesn't compress based on the file extension — this needs to be fixed.
- Biological note: `-f 4` keeps reads where that specific read is unmapped; it needs to be decided whether only pairs where both mates are non-host should be kept.

---

## `workflow/rules/meta_genome.smk`

### Settings

- Gets the list of samples.
- Selects the classification method and Kraken database.
- Copies paths/prefix.
- Restricts the level wildcard to `S` only. So with the file as it stands, genus/phylum outputs aren't possible, even though the comment mentions them.
- A Python function that returns the input path for a given wildcard. It expects `.fastq`, whereas the previous rule declares `.fastq.gz` and a different output directory.

### Kraken2 branch

- This block is selected if the config specifies `kraken2`.
- Stops the workflow early if there's no database path.
- Defines `kraken2_classify`.
- Calls the `INPUT_R1` function, passing it a specific `sra_id`.
- Marks the report and the raw classification output as `temp()`. Snakemake can delete these files after downstream rules run, which conflicts with wanting to keep the reports as results — `temp()` should be removed for archival outputs.
- Runs Kraken2 on a single file, so this is single-end syntax.

### Bracken

- Requires the Kraken report — this is the Kraken → Bracken dependency.
- Produces the abundance report per sample and rank.
- Uses `READ_LEN`, but that variable is only defined in the root `Snakefile`; if the `.smk` file is run on its own, it won't exist. Using `config["read_length"]` directly would be more robust.
- Sets an abundance threshold of 10 reads, but this parameter isn't in the YAML and isn't documented in the README.
- Runs Bracken with the database, input report, output, read length, rank, and threshold.

### Merging

- Builds the list of all species reports.
- Produces a single combined `gallus_gallus_report.tsv`.
- Calls `kreport2mpa.py`. It should be checked whether this utility actually accepts Bracken output in this format — it's more commonly applied to a Kraken report.

### MetaPhlAn alternative

- This branch is only selected when `algorithm: metaphlan`.
- Uses `INPUT_R2`, but no such function exists in the file.
- The whole branch produces a single shared output for all samples and has no wildcard, so it doesn't actually implement multi-sample profiling.
- This is a draft and shouldn't be run without further work.

---

## `workflow/rules/transcriptome.smk`

### Paths

- Builds the fasta path as `ref_dir/{taxid}.fa`. By its name this is the transcriptome fasta, but the directory `ref_dir` is called `genome` — it should be confirmed that it actually holds transcript sequences and not genomic fasta.
- Again fetches the list of SRA IDs, even though the variable is named in the singular.
- Uses `{SRA_ID}` in uppercase and `input_dir`, which in the config points to the non-host metagenome reads. Host RNA-seq needs a separate `host_transcriptome_input_dir` and the `{sra_id}` wildcard.

### `kallisto_index`

- The transcript fasta.
- The binary Kallisto index.
- Hard-codes k=31 instead of using `config["kallisto"]["kmer_size"]`.
- Calls `kallisto index`.

### `kallisto_quant`

- The index and the two mate read files.
- The output directory, which will contain `abundance.tsv` and Kallisto's auxiliary files.
- `-b 30` bootstraps.
- Checks whether R2 exists. But since R2 is declared as a required input, Snakemake won't even start the job if it's missing — so the single-end branch is effectively unreachable.
- The paired-end quantification.
- The intended single-end quantification.
- This file is missing `import os`, so the R2-existence check will fail at interpretation time.

---

## `workflow/rules/meta_transcriptome.smk`

This is a separate draft that isn't currently included in the root `Snakefile`.

### QC rule

- Looks for `config.yaml` in the current directory, not `config/config.yaml`.
- Expects `input_fastq`, which doesn't exist in the current config.
- Produces a single shared output with no sample wildcard.
- Runs single-end fastp with Q30 and a minimum length of 50.

### Supposed MetaPhlAn step

- The input and output names look reasonable.
- However, it runs `fastp` again and references `input.fastq`, `output.processed_fastq`, `output.html_report`, none of which exist in this rule.
- In other words, metatranscriptomic taxonomic profiling **isn't implemented**.

### HUMAnN

- Describes the intended input/output of the functional profile.
- Contains a second, redundant `shell:` — a syntax error.
- The intended HUMAnN command: a functional abundance table.

### HUMAnN postprocessing

- `rename_human` renames pathway identifiers to UniRef90 names.
- `regroup_table_humman` regroups features by MetaCyc. The name `humman` is a historical typo.
- `final_humann_table` again uses `humann_rename_table`, but passes `--groups metacyc`, which doesn't match the purpose of the rename command — the CLI and the intended goal of this step should be checked.

### Cleanup and target

- `clean_temporary_files` removes the processed fastq and creates a marker.
- Its output isn't required anywhere, so the rule won't necessarily run.
- `rule all` requests the MetaPhlAn and final HUMAnN outputs, but as things stand, the upstream MetaPhlAn/HUMAnN rules aren't valid.