# Contributing

This project and everyone participating in it is governed by the [CONTRIBUTING.md](./master/docs/CONTRIBUTING.md).
Follow existing style, formatting, and naming conventions for the file you are modifying and for the entire project.
Please review the following guidelines before you get started


## Reporting bugs and proposing features

- Search [existing issues](../../issues) before opening a new one.
- For a bug report, include: what you ran, what you expected, what happened instead, and the relevant log file (pipeline logs are written under `{output_dir}/logs/`).
- For a feature or pipeline-stage proposal, describe the input/output of the new stage and which existing rule it should connect to.

## Code Style

- Follow the PEP 8 standards for Python
- Variable and function names should be unique and reflect their purpose
- Class names should start with a capital letter
- Add comments to complex code.
- Use flake8 and black to check and automatically format code

##  Function Annotation

- The Python typing module is used to annotate types for functions and variables.
- Annotate the types of all functions, methods, arguments, and return values, as well as variables, where possible and appropriate.

- Using:
  - **str**, **int**, **float**, and **bool** to annotate simple data types.
  - **List**, **Tuple**, **Dict**, and **Set** from the typing module to annotate collections.
  - Union to indicate that a variable can have one of several types.
  - Docstrings to describe the purpose of a function, its arguments, and its return value.


## ArgParse annotation

- Use ArgParse for all scripts that need the user to interact with the command line.
- Using:
  - Create an ArgumentParser object.
  - Describe what the script does in the description
  - Use the *add_argument()* method to define the arguments that your script will accept.

    The options for *add_argument()* :
    - *name or flags*: The name of the argument or a list of flags (e.g., -o, --output).
    - *help*: A short description of the argument to be displayed in the help.
    - *type*: The expected data type of the argument (e.g., str, int, float).
    - *default*: The default value of the argument.
    - *required*: Specify True if the argument is required.
    - *choices*: A list of acceptable values for the argument.
    - *action*: Specifies the action to perform when this argument is detected on the command line.
      
  - Calling the parse_args() method to parse command line arguments

## Base Packages

The following basic packages are required to work on the project: 

- python (version 3.9 or higher)
- pytest (for running tests)
- flake8 (for code style checking)
- black (for automatic code formatting)
- requests (for working with HTTP requests)
- pandas (for data analysis)
- numpy (for data analysis)

The minimum list of packages and their versions that are required for the project can be found in requirements/minimal.txt

## Installing packages

- Every package is installed conda, save any installation code you add to the repository (`requirements/installations.sh`).
- Conda dependencies and their versions are listed in `enviromental.yaml` at the repository root. 

## Snakemake conventions

- Each pipeline stage are in its own rule file under `workflow/rules/` and is pulled into the root `Snakefile` via `include:`.
- Give every rule a `log:` directive pointing into `{output_dir}/logs/`.

## Folder structure

This is the actual structure of the repository — keep new code inside it rather than introducing parallel folders:

- `workflow/rules/`: Snakemake rule files, one per pipeline stage (preprocessing, host/metagenome split, metagenome classification, transcriptome quantification, metatranscriptome, genome download).
- `scripts/`:  Python/Bash scripts invoked by the pipeline or run manually.
- `CCA_analysis/`: R scripts for the CCA and downstream statistical analysis.
- `config/`: `config.yaml` and any other run configuration.
- `data/`: input data references (SRA ID lists, reference genomes), large/raw data files are excluded via `.gitignore`.
- `tests/`: tests for `scripts/` and for the Snakemake rules. Uses pytest, test filenames start with `test_`.
- `docs/`: documentation, including this file.
- `requirements/`: minimal package lists (`minimal.txt`, `tests.txt`) and install scripts.
- Code generated during development, testing, or builds must be in `.gitignore`, not committed.

## Running tests locally

```
pip install -r requirements/tests.txt
pytest tests/
```

## Pull requests

- Describe the changes made in the PR description clearly and in detail.
- Link the issue(s) the PR addresses.
- Verify that all tests pass, and that `snakemake -n` (dry-run) succeeds if you touched any rule file, before requesting review.
- Update documentation (README, `docs/PIPELINE.md`, this file) if your change affects usage, configuration, or folder structure.