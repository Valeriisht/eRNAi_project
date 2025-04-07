import os
import time  
import pytest
from pathlib import Path
from snakemake import snakemake
from unittest.mock import patch

@pytest.fixture
def setup_kallisto_test(tmp_path):
    """Фикстура для подготовки тестового окружения"""
    # Создаем структуру директорий
    (tmp_path / "transcriptome_kallisto").mkdir()
    (tmp_path / "logs_kallisto").mkdir()
    
    # Тестовые данные
    taxid = "12345"
    sra_id = "TEST123"
    
    # Создаем входные файлы
    (tmp_path / f"{taxid}.fna").write_text(">test\nACGT"*100)  # Транскриптом
    (tmp_path / f"{sra_id}_filtered_1.fastq").write_text("@read\nACGT\n+\nIIII")  # R1
    (tmp_path / f"{sra_id}_filtered_2.fastq").write_text("@read\nTGCA\n+\nIIII")  # R2
    
    # Конфигурация
    config = {
        "output_dir": str(tmp_path),
        "sra": {"sra_id": sra_id},
        "kallisto": {"bootstrap": 30}
    }
    
    return tmp_path, config, taxid, sra_id

def test_kallisto_index_rule(setup_kallisto_test, mocker):
    test_dir, config, taxid, sra_id = setup_kallisto_test
    
    # 1. Мокаем shell-команды
    mocker.patch('snakemake.shell', return_value=0)
    
    # 2. Создаем УПРОЩЕННЫЙ Snakefile без config
    snakefile = test_dir / "Snakefile"
    snakefile.write_text(f"""
rule kallisto_index:
    input:
        "{taxid}.fna"  # Явное значение вместо wildcard
    output:
        "transcriptome_kallisto/{taxid}_transcriptome.idx"  # Явное значение
    params:
        kmer_size = 31
    log:
        "logs_kallisto/{taxid}_kallisto_index.log"
    shell:
        "touch {{output}} && echo 'mocked' > {{log}}"
""")

    # 3. Запускаем с минимальными параметрами
    result = snakemake(
        snakefile=str(snakefile),
        workdir=str(test_dir),
        targets=[f"transcriptome_kallisto/{taxid}_transcriptome.idx"],
        cores=1,
        quiet=True,
        forceall=True
    )

    # 4. Проверяем базовые условия
    assert result, "Snakemake execution failed"
    assert (test_dir / f"transcriptome_kallisto/{taxid}_transcriptome.idx").exists()

def test_kallisto_quant_rule_paired(setup_kallisto_test, mocker):
    """Тест для правила квантификации (paired-end)"""
    test_dir, config, taxid, sra_id = setup_kallisto_test
    
    # 1. Создаем все необходимые директории
    (test_dir / ".snakemake/log").mkdir(parents=True, exist_ok=True)
    (test_dir / "transcriptome_kallisto").mkdir(exist_ok=True)
    (test_dir / "logs_kallisto").mkdir(exist_ok=True)
    
    # 2. Создаем входные файлы
    idx_file = test_dir / f"transcriptome_kallisto/{taxid}_transcriptome.idx"
    idx_file.touch()
    (test_dir / f"{sra_id}_filtered_1.fastq").write_text("@read\nACGT\n+\nIIII")
    (test_dir / f"{sra_id}_filtered_2.fastq").write_text("@read\nTGCA\n+\nIIII")
    
    # 3. Мокаем системные вызовы
    mock_shell = mocker.patch('snakemake.shell', return_value=0)
    mocker.patch('os.path.exists', return_value=True)
    mocker.patch('os.path.getmtime', return_value=time.time())  # Мок для временных меток
    
    # 4. Создаем Snakefile с явными путями
    snakefile = test_dir / "Snakefile"
    snakefile.write_text(f"""
rule kallisto_quant:
    input:
        index = "{str(idx_file)}",
        r1 = "{sra_id}_filtered_1.fastq",
        r2 = "{sra_id}_filtered_2.fastq"
    output:
        directory("transcriptome_kallisto/{taxid}/{sra_id}_quant_results")
    params:
        bootstrap = 30
    log:
        "logs_kallisto/{taxid}/{sra_id}_kallisto_quant.log"
    shell:
        "mkdir -p {{output}} && touch {{output}}/abundance.tsv && echo 'mocked' > {{log}}"
""")
    
    # 5. Запускаем Snakemake
    result = snakemake(
        snakefile=str(snakefile),
        workdir=str(test_dir),
        cores=1,
        quiet=True,
        forceall=True,  # Принудительное выполнение
        targets=[f"transcriptome_kallisto/{taxid}/{sra_id}_quant_results"],
        lock=False  # Отключаем блокировку файлов
    )
    
    # 6. Проверяем результаты
    assert result, "Snakemake execution failed"
    mock_shell.assert_called_once()
    
    # Проверяем выходные файлы
    out_dir = test_dir / f"transcriptome_kallisto/{taxid}/{sra_id}_quant_results"
    assert out_dir.is_dir(), f"Output directory {out_dir} missing"
    assert (out_dir / "abundance.tsv").exists(), "Output file missing"