import pytest
from snakemake import snakemake

@pytest.fixture(scope="session")
def snakemake_runner():
    def run(rule, config=None):
        return snakemake(
            snakefile="Snakefile",
            targets=[rule],
            config=config or {},
            workdir=".",
            quiet=True,
            cores=1
        )
    return run