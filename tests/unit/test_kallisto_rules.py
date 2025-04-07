import os
import pytest
from pathlib import Path
from snakemake import snakemake
import tempfile

@pytest.fixture
def setup_kallisto_test(tmp_path):
    """Фикстура для тестов Kallisto"""
    # test structure
    test_dir = tmp_path / "kallisto_test"
    test_dir.mkdir()
    
    # config file
    config = {
        "output_dir": str(test_dir),
        "sra": {"sra_id": "TEST123"},
        "kallisto": {"bootstrap": 30}
    }
    
    # tests file
    taxid = "12345"
    (test_dir / f"{taxid}.fna").write_text(">test\nACGT"*100)
    (test_dir / f"TEST123_filtered_1.fastq").write_text("@read\nACGT\n+\nIIII")
    (test_dir / f"TEST123_filtered_2.fastq").write_text("@read\nTGCA\n+\nIIII")
    
    return test_dir, config, taxid

def test_kallisto_index_rule(setup_kallisto_test, mocker):
    """Transcriptome index rule"""
    test_dir, config, taxid = setup_kallisto_test
    mock_shell = mocker.patch('snakemake.shell')
    
    snakefile = test_dir / "Snakefile"
    snakefile.write_text(f"""
configfile: "config.yaml"

rule kallisto_index:
    input:
        transcriptome = config["output_dir"] + "/{{taxid}}.fna"
    output:
        transcriptome_index = config["output_dir"] + "/transcriptome_kallisto/{{taxid}}_transcriptome.idx"
    params:
        kmer_size = 31
    log:
        "logs_kallisto/{{taxid}}_kallisto_index.log"
    shell:
        "kallisto index -i {{output.transcriptome_index}} -k {{params.kmer_size}} {{input.transcriptome}} > {{log}} 2>&1"
    """)
    
    # Create config.yaml
    (test_dir / "config.yaml").write_text(f"""
output_dir: {config['output_dir']}
sra:
  sra_id: {config['sra']['sra_id']}
kallisto:
  bootstrap: {config['kallisto']['bootstrap']}
    """)
    
    # Launch 
    result = snakemake(
        snakefile=str(snakefile),
        workdir=str(test_dir),
        configfiles=[str(test_dir / "config.yaml")],
        targets=[f"transcriptome_kallisto/{taxid}_transcriptome.idx"],
        cores=1
    )
    
    assert result
    mock_shell.assert_called_once()
    call_args = mock_shell.call_args[0][0]
    assert f"-i {test_dir}/transcriptome_kallisto/{taxid}_transcriptome.idx" in call_args
    assert "-k 31" in call_args
    assert f">{test_dir}/logs_kallisto/{taxid}_kallisto_index.log" in call_args

def test_kallisto_quant_rule_paired(setup_kallisto_test, mocker):
    """Тест правила квантификации (paired-end)"""
    test_dir, config, taxid = setup_kallisto_test
    mock_shell = mocker.patch('snakemake.shell')
    
    snakefile = test_dir / "Snakefile"
    snakefile.write_text(f"""
configfile: "config.yaml"

rule kallisto_quant:
    input:
        index = config["output_dir"] + "/transcriptome_kallisto/{{taxid}}_transcriptome.idx",
        r1 = config["output_dir"] + "/{{SRA_ID}}_filtered_1.fastq",
        r2 = config["output_dir"] + "/{{SRA_ID}}_filtered_2.fastq"
    output:
        directory(config["output_dir"] + "/transcriptome_kallisto/{{taxid}}/{{SRA_ID}}_quant_results")
    params:
        bootstrap = config["kallisto"]["bootstrap"]
    log:
        "logs_kallisto/{{taxid}}/{{SRA_ID}}_kallisto_quant.log"
    run:
        if os.path.exists(input.r2):
            shell(
                "kallisto quant -i {{input.index}} -o {{output}} -b {{params.bootstrap}} {{input.r1}} {{input.r2}} > {{log}} 2>&1"
            )
    """)
    
    # Create index
    (test_dir / "transcriptome_kallisto").mkdir()
    (test_dir / f"transcriptome_kallisto/{taxid}_transcriptome.idx").touch()
    
    result = snakemake(
        snakefile=str(snakefile),
        workdir=str(test_dir),
        configfiles=[str(test_dir / "config.yaml")],
        targets=[f"transcriptome_kallisto/{taxid}/TEST123_quant_results"],
        cores=1
    )
    
    assert result
    mock_shell.assert_called_once()
    call_args = mock_shell.call_args[0][0]
    assert f"-i {test_dir}/transcriptome_kallisto/{taxid}_transcriptome.idx" in call_args
    assert f"-o {test_dir}/transcriptome_kallisto/{taxid}/TEST123_quant_results" in call_args
    assert f"-b {config['kallisto']['bootstrap']}" in call_args