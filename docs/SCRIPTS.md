# Review of helper scripts, environments, and tests

## `scripts/download_metagenome.py`

This Python script downloads **reference bacterial genomes**. It's needed after CCA, for the taxa identified in the ordination, to obtain FASTA genomes for k-mer/sequence analysis.

A few notes on the code follow.

### Imports and logging

- `gzip` opens `.gz` FASTA files without manual decompression.
- `tempfile` creates an automatically cleaned-up temporary directory.
- `subprocess` runs the external `ncbi-genome-download` command.
- `dataclass` creates a compact configuration class.

### `GenomeDownloadConfig`

`@dataclass` automatically generates the class constructor. Fields:

| Field | Default value |
| --- | --- |
| `assembly_level` | complete assemblies only |
| `refseq_category` | RefSeq reference genomes only |
| `file_format` | FASTA |
| `parallel` | 3 parallel downloads |
| `retries` | 3 attempts |

A parameter template passed later into the function.

### `check_existing_genome`

- Arguments: genus name and output folder.
- Builds the expected filename, e.g. `output_downloaded_metagenome/Escherichia_reference.fna.gz`.
- Returns a `Path` object if the file exists, otherwise `None`.

The script stores **one genome per genus**, not one genome per species.

Issue: for CCA taxa at species level this is a very rough simplification — different `Bacteroides` species, for example, could end up represented by the wrong genome.

### `verify_genome_file`

- Opens the gzip FASTA in text mode.
- Reads only the first header line.
- Considers the genome a match if the genus name appears in the header.
- Logs an error otherwise.

Issue: this doesn't guarantee taxonomic correctness — the header may use an abbreviated name, or the genus may appear in a different field.

### `download_genomes`

Arguments: a collection/list of genera, an output directory, and a configuration.

- Creates the target directory.
- Builds the result as a `{genus: downloaded_path}` mapping.
- Iterates over the taxa: reuses a previously saved valid file, or adds the genus to `genera_to_download`.
- Returns early if everything is already downloaded.
- Creates a temporary directory that is removed at the end of the block.
- Builds the external program's arguments as an argument list rather than a shell string — safer than constructing a shell command.
- Runs `ncbi-genome-download` with `check=True`, which raises an exception on a nonzero exit code.
- Expects a specific output structure from that program: iterates over `*_genomic.fna.gz`, checks the header, and copies the first matching file into the permanent directory.
- Catches errors and only logs them rather than re-raising — so the caller must check for an incomplete `results`.
- Returns the dictionary of successfully found genomes.

---

## `scripts/run_download_metagenome.py`

A short driver for the module above.

- Imports `download_genomes`.
- Computes a directory two levels above the current working directory.
  Issue: the result depends on where the script is run from.
- Expects `species.csv` in that computed directory.
- Reads the CSV, treats the first column as the index, strips whitespace after delimiters, and expects `Genus`, `CCA1`, `CCA2` columns.
- Sets the output folder.
- Strips leading/trailing whitespace from genus names.
- Converts the `Genus` column into a Python list.
- Downloads the genomes.

`species.csv` contains CCA coordinates; only the genus names are actually used here.

---

## `scripts/find_shared_kmers.sh`

This script assumes DSK has already produced `kmer-count`-style tables, and finds **exact shared k-mers** between the host genome and each metagenome file.

### Setup

- Selects the Bash interpreter.
- Sets the expected directory of input metagenome DSK tables.
- Sets the DSK table for the chicken genome.
- Sets the minimum occurrence count for a k-mer in the host genome to 2.

### `process_file`

- Takes a single metagenome DSK table as positional argument `$1`.
- Extracts the sample name from the basename, stripping only the `.txt` suffix.
- Builds the final output filename plus two temporary filenames.
- If the cleaned host table doesn't exist yet, uses `awk '{print $1}'` to take the first column (the k-mers themselves) from the host DSK table; `||` means "otherwise."
- Similarly extracts k-mers from the current metagenome table.
- Uses `grep -Fxf`: `-F` for literal strings, `-x` for full-line match, `-f` for patterns from the host table — the result is the exact intersection of k-mers.
- Checks whether a frequency filter is needed.
- The first `awk` reads the host k-mer/count table into an associative array; the second input is the intersection lines from the metagenome. It prints the k-mer and the **host genome count** if that count is ≥ 2. The metagenome count is neither filtered nor saved to the final file.
- Counts the final rows, writes a summary, and returns success.
- On failure, writes an error log and returns failure.

### Loop

- Exports the function/variables for GNU Parallel, but `parallel` is never actually called — a leftover from a possible parallel version.
- Iterates only over `*_k15.h5.txt` filenames: k=15 is hardcoded into the filename.
- Retries processing up to three times, waiting 10 seconds between attempts.
- Removes temporary files if the last command in the loop succeeded. If an earlier file failed, this check may not reflect that correctly — an explicit error counter or `set -euo pipefail` would be safer.

### Interpreting the result

`*_shared_kmers.txt` answers the question: "which short DNA strings occur in both the host genome and the selected metagenome table?"

---

### `requirements/installations.sh`

Each line is a separate interactive `conda install` command.
This file should be replaced with a single lockable `environment.yml` with pinned versions. It can be removed.

### `requirements/minimal.txt`

A minimal Python development list: Python, pytest, flake8, black, requests, pandas, numpy. It doesn't install the conda bioinformatics programs and isn't sufficient for the full pipeline.

### `requirements/tests.txt`

A separate conda YAML: `snakemake-tests` environment, Python 3.10, pytest, Snakemake, Kallisto, pytest-mock.

---

## Tests: to do

### `tests/__init__.py`

Makes `tests/` a package.

### `tests/conftest.py`

- Provides the `mocker` fixture via the pytest-mock plugin.

### `tests/unit/conftest.py`



### `tests/unit/test_kallisto_rules.py`



### `tests/unit/test_sra_rules.py`

- These are useful smoke tests documenting the expected commands, but not integration tests. They won't catch syntax errors, missing paths, or incorrect connections between the real Snakemake rules.