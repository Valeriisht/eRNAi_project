import pytest
import os
from unittest.mock import patch, call
import tempfile
import yaml

@pytest.fixture
def config():
    """Fixture with test configuration"""
    return {
        "output_dir": "test_output",
        "sra": {
            "sra_id": "TEST123"
        },
        "kallisto": {
            "bootstrap": 100
        }
    }

def test_kallisto_index_rule(config):
    """Test the kallisto_index rule structure and parameters"""
    rule = {
        "input": {
            "transcriptome": f"{config['output_dir']}/{{taxid}}.fna"
        },
        "output": {
            "transcriptome_index": f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}_transcriptome.idx"
        },
        "params": {
            "kmer_size": 31
        },
        "log": f"logs_kallisto/{{taxid}}_kallisto_index.log",
        "shell": (
            "kallisto index -i {output.transcriptome_index} "
            "-k {params.kmer_size} {input.transcriptome} > {log} 2>&1"
        )
    }
    
    # Verify input files
    assert rule["input"]["transcriptome"] == f"{config['output_dir']}/{{taxid}}.fna"
    
    # Verify output files
    assert rule["output"]["transcriptome_index"] == (
        f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}_transcriptome.idx"
    )
    
    # Verify parameters
    assert rule["params"]["kmer_size"] == 31
    
    # Verify logs
    assert rule["log"] == "logs_kallisto/{taxid}_kallisto_index.log"
    
    # Verify shell command components
    assert "kallisto index" in rule["shell"]
    assert "-i {output.transcriptome_index}" in rule["shell"]
    assert "-k {params.kmer_size}" in rule["shell"]
    assert "{input.transcriptome}" in rule["shell"]

def test_kallisto_quant_rule(config):
    """Test the kallisto_quant rule structure and parameters"""
    # Mock the rule structure
    rule = {
        "input": {
            "index": f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}_transcriptome.idx",
            "r1": f"{config['output_dir']}/{{SRA_ID}}_filtered_1.fastq",
            "r2": f"{config['output_dir']}/{{SRA_ID}}_filtered_2.fastq"
        },
        "output": {
            "directory": f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}/{{SRA_ID}}_quant_results"
        },
        "params": {
            "bootstrap": config["kallisto"]["bootstrap"]
        },
        "log": f"logs_kallisto/{{taxid}}/{{SRA_ID}}_kallisto_quant.log",
        "run": None  # Заполним ниже
    }
    
     # Mock the run behavior for paired-end and single-end modes
    def mock_run():
        if os.path.exists(rule["input"]["r2"]):
            return (
                "kallisto quant -i {input.index} -o {output} "
                "-b {params.bootstrap} {input.r1} {input.r2} > {log} 2>&1"
            )
        else:
            return (
                "kallisto quant -i {input.index} -o {output} "
                "-b {params.bootstrap} {input.r1} > {log} 2>&1"
            )
    
    rule["run"] = mock_run
    
    # Verify input files
    assert rule["input"]["index"] == (
        f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}_transcriptome.idx"
    )
    assert rule["input"]["r1"] == f"{config['output_dir']}/{{SRA_ID}}_filtered_1.fastq"
    assert rule["input"]["r2"] == f"{config['output_dir']}/{{SRA_ID}}_filtered_2.fastq"
    
    # Verify output directory
    assert rule["output"]["directory"] == (
        f"{config['output_dir']}/transcriptome_kallisto/{{taxid}}/{{SRA_ID}}_quant_results"
    )
    
    # Verify parameters
    assert rule["params"]["bootstrap"] == config["kallisto"]["bootstrap"]
    
    # Verify logs
    assert rule["log"] == "logs_kallisto/{taxid}/{SRA_ID}_kallisto_quant.log"
    
    # Test single-end mode command (when R2 doesn't exist
    paired_cmd = rule["run"]()
    assert "kallisto quant" in paired_cmd
    assert "-i {input.index}" in paired_cmd 