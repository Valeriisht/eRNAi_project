import pytest
from pathlib import Path


@pytest.fixture
def setup_kallisto_test(tmp_path):
    test_dir = tmp_path / "kallisto_test"
    test_dir.mkdir()
    

    (test_dir / "12345.fna").write_text(">test\nACGT"*100)
    (test_dir / "TEST123_filtered_1.fastq").write_text("@read\nACGT\n+\nIIII")
    (test_dir / "TEST123_filtered_2.fastq").write_text("@read\nTGCA\n+\nIIII")
    
    # config
    config = {
        "output_dir": str(test_dir),
        "sra": {"sra_id": "TEST123"},
        "kallisto": {"bootstrap": 30}
    }
    
    return test_dir, config, "12345"


@pytest.fixture
def mock_exists(mocker):
    return mocker.patch('os.path.exists')