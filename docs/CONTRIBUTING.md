# Contributing

This project and everyone participating in it is governed by the [CONTRIBUTING.md](./master/docs/CONTRIBUTING.md).
Follow existing style, formatting, and naming conventions for the file you are modifying and for the entire project.
Please review the following guidelines before you get started

## Code Style

- Follow the PEP 8 standards for Python
- Variable and function names should be unique and reflect their purpose
- Class names should start with a capital letter
- Add comments to complex code.
- Use flake8 and black to check and automatically format code

##  Function Annotation

- The Python typing module is used to annotate types for functions and variables.
- Use the argparse module to process command line arguments

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

## Install packeges 

- Each package is installed using conda. Save the installation code in repository. 

- The repository contains a separate YAML files that lists the dependencies and the versions used. A list of Conda packages and channels to install is contained in the environment.yaml file.

## SnakeMake wraps

- The enviromental.yaml file for each rule must be stacked in the Snakefile so that Snakemake can create an isolated conda environment with the necessary dependencies to run that rule.

## Folder structure

Use the following folder structure

- src/: Framework source code:
  - Contains modules and packages that implement the functionality of the framework. Each module should have a unique name and contain docstrings.

- tests/: Tests.
  - Contains tests for every module in src/.
  - The pytest library is used to write and execute the tests. Test filenames start with test_.

- docs/: Documentation
  - Includes a file with the contributor rules 

- examples/: Usage examples.
  - Contains code examples that demonstrate how to use the framework.

- requirements/: packages 
  - contains a set of minimum required packages 

- scripts/: Basic project scripts for execution.

- Code generated during development, testing, build/compile should be in .gitignore

## Execute through a Python or Bash layer

- Use Python or Bash scripts to execute the framework

## Pull Requests

- Describe changes made in the PR description clearly and in detail.
- Verify that all tests pass before submitting the PR.
- Update documentation if you make changes to the code.
  
For function annotation, use the Callable object from the typing module.



