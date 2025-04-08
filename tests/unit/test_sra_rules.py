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
            "sra_id": "TEST123",
            "thread": 4,
            "paired": True
        },
        "fastp": {
            "thread": 2,
            "qualified_quality_phred": 20,
            "min_length": 50,
            "detect_adapters": True
        }
    }

@pytest.fixture
def mock_os():
    with patch('os.makedirs') as mock_makedirs, \
         patch('os.path.join') as mock_join:
        yield mock_makedirs, mock_join

def test_prefetch_data_rule(config, mock_os):
    mock_makedirs, mock_join = mock_os
    
    # Mock rule structure
    rule = {
        "output": {
            "sra_file": "results/sra/{SRA_ID}.sra"
        },
        "params": {
            "sra_id": "{SRA_ID}"
        },
        "log": config["output_dir"] + "/logs/{SRA_ID}_prefetch.log",
        "shell": "prefetch {params.sra_id} --output-file {output.sra_file} > {log} 2>&1"
    }
    
    # Verify output parameters
    assert "results/sra/{SRA_ID}.sra" in rule["output"]["sra_file"]
    
    # Verify parameters
    assert rule["params"]["sra_id"] == "{SRA_ID}"
    
    # Verify logs
    assert f"{config['output_dir']}/logs/{{SRA_ID}}_prefetch.log" in rule["log"]
    
    # Check shell
    assert "prefetch {params.sra_id} --output-file {output.sra_file}" in rule["shell"]

def test_download_data_rule(config):
    # Mock rules
    rule = {
        "input": {
            "sra_file": "results/sra/{SRA_ID}.sra"
        },
        "output": {
            "f1": f"{config['output_dir']}/raw_{{SRA_ID}}_1.fastq",
            "r1": f"{config['output_dir']}/raw_{{SRA_ID}}_2.fastq"
        },
        "params": {
            "sra_id": "{SRA_ID}",
            "threads": config["sra"]["thread"],
            "paired": config["sra"]["paired"]
        },
        "log": f"{config['output_dir']}/logs/{{SRA_ID}}_download.log",
        "shell": """
        set -euo pipefail
        strace -e trace=file fasterq-dump {input.sra_file} \
        --outdir {OUTPUT_DIR} \
        --split-files  \
        --threads {params.threads} > {log} 2>&1

        # Переименование файлов
        mv {OUTPUT_DIR}/{params.sra_id}_1.fastq {output.f1}
        mv {OUTPUT_DIR}/{params.sra_id}_2.fastq {output.r1}
        """
    }
    
    # Verify input files
    assert rule["input"]["sra_file"] == "results/sra/{SRA_ID}.sra"
    
    # Verify output files
    assert rule["output"]["f1"] == f"{config['output_dir']}/raw_{{SRA_ID}}_1.fastq"
    assert rule["output"]["r1"] == f"{config['output_dir']}/raw_{{SRA_ID}}_2.fastq"
    
    # Verify parameters
    assert rule["params"]["sra_id"] == "{SRA_ID}"
    assert rule["params"]["threads"] == config["sra"]["thread"]
    assert rule["params"]["paired"] == config["sra"]["paired"]
    
    # Check shell
    assert "fasterq-dump" in rule["shell"]
    assert "--split-files" in rule["shell"]
    assert "--threads {params.threads}" in rule["shell"]
    assert "mv {OUTPUT_DIR}/{params.sra_id}_1.fastq {output.f1}" in rule["shell"]

def test_process_paired_data_rule(config):
    # Mock rules
    rule = {
        "input": {
            "f1": f"{config['output_dir']}/raw_{{SRA_ID}}_1.fastq",
            "r1": f"{config['output_dir']}/raw_{{SRA_ID}}_2.fastq"
        },
        "output": {
            "filtered_f1": f"{config['output_dir']}/{{SRA_ID}}_filtered_1.fastq",
            "filtered_r2": f"{config['output_dir']}/{{SRA_ID}}_filtered_2.fastq",
            "report_json": f"{config['output_dir']}/{{SRA_ID}}_fastp_report.json"
        },
        "params": {
            "threads": config["fastp"]["thread"],
            "quality_threshold": config["fastp"]["qualified_quality_phred"],
            "min_length": config["fastp"]["min_length"],
            "detect_adapters": config["fastp"]["detect_adapters"]
        },
        "log": f"{config['output_dir']}/logs/{{SRA_ID}}_fastp_paired.log",
        "shell": """
        fastp -i {input.f1} -I {input.r1} \
        -o {output.filtered_f1} -O {output.filtered_r2} \
        --thread {params.threads} \
        --detect_adapter_for_pe \
        -q {params.quality_threshold} \
        --length_required {params.min_length} \
        --json {output.report_json} > {log} 2>&1

        rm {input.f1} {input.r1}
        """
    }
    
    # Verify input files
    assert rule["input"]["f1"] == f"{config['output_dir']}/raw_{{SRA_ID}}_1.fastq"
    assert rule["input"]["r1"] == f"{config['output_dir']}/raw_{{SRA_ID}}_2.fastq"
    
    # Verify output files
    assert rule["output"]["filtered_f1"] == f"{config['output_dir']}/{{SRA_ID}}_filtered_1.fastq"
    assert rule["output"]["filtered_r2"] == f"{config['output_dir']}/{{SRA_ID}}_filtered_2.fastq"
    assert rule["output"]["report_json"] == f"{config['output_dir']}/{{SRA_ID}}_fastp_report.json"
    
    # Verify parameters
    assert rule["params"]["threads"] == config["fastp"]["thread"]
    assert rule["params"]["quality_threshold"] == config["fastp"]["qualified_quality_phred"]
    assert rule["params"]["min_length"] == config["fastp"]["min_length"]
    assert rule["params"]["detect_adapters"] == config["fastp"]["detect_adapters"]
    
    # Command shell
    assert "fastp -i {input.f1} -I {input.r1}" in rule["shell"]
    assert "-o {output.filtered_f1} -O {output.filtered_r2}" in rule["shell"]
    assert "--thread {params.threads}" in rule["shell"]
    assert "-q {params.quality_threshold}" in rule["shell"]
    assert "--length_required {params.min_length}" in rule["shell"]
    assert "--json {output.report_json}" in rule["shell"]
    assert "rm {input.f1} {input.r1}" in rule["shell"]

def test_directory_creation(config, mock_os):
    mock_makedirs, mock_join = mock_os
    
    # Mock the rule structure
    os.makedirs(config["output_dir"], exist_ok=True)
    os.makedirs(os.path.join(config["output_dir"], "logs"), exist_ok=True)
    
    # Check mock calls
    mock_makedirs.assert_has_calls([
        call(config["output_dir"], exist_ok=True),
        call(os.path.join(config["output_dir"], "logs"), exist_ok=True)
    ])

@pytest.fixture
def temp_config_file(config):
    # Create tmp config
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        yaml.dump(config, f)
        return f.name

def test_config_loading(temp_config_file):
    # Test confif load
    with open(temp_config_file) as f:
        loaded_config = yaml.safe_load(f)
    
    assert loaded_config["output_dir"] == "test_output"
    assert loaded_config["sra"]["sra_id"] == "TEST123"
    assert loaded_config["fastp"]["thread"] == 2
    
    # Delete tmp file
    os.unlink(temp_config_file)