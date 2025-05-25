import gzip
import tempfile
from pathlib import Path
from typing import Dict, Set, Optional
import subprocess
from dataclasses import dataclass
import logging 


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class GenomeDownloadConfig:
    """Config for genome download
    
    Attributes:
        assembly_level: Genome assembly level ('complete', 'chromosome', etc.)
        refseq_category: RefSeq category ('reference', 'representative')
        file_format: Download format ('fasta', 'genbank', etc.)
        parallel: Number of parallel downloads
        retries: Number of retry attempts
    
    """
    assembly_level: str = "complete"
    refseq_category: str = "reference"
    file_format: str = "fasta"
    parallel: int = 3
    retries: int = 3


def check_existing_genome(genus: str, genome_dir: Path) -> Optional[Path]:

    """
    Check if genome file already exists for the specified genus.
    
    Args:
        genus: Genus name (e.g., 'Escherichia')
        genome_dir: Directory containing genome files
        
    Returns:
        Path to existing genome file if found, None otherwise
    
    """
    genome_path = genome_dir / f"{genus}_reference.fna.gz"
    return genome_path if genome_path.exists() else None 


def verify_genome_file(genome_file: Path, genus: str) -> bool:
    """Verify that the genome file actually contains DNA sequences for the expected genus.
    
    Args:
        genome_file: Path to the genome file
        genus: Expected genus name
        
    Returns:
        True if verification passes, False otherwise
    """

    try:
        with gzip.open(genome_file, "rt") as f:
            header = f.readline()
            return genus.lower() in header.lower()
    
    except Exception as e:
        logger.error(f"Error verifying {genome_file}: {e}")


def download_genomes(
    genera: Set[str],
    output_dir: Path,
    config: GenomeDownloadConfig = GenomeDownloadConfig()
) -> Dict[str, Path]:
    """Download reference genomes from NCBI for specified bacterial genera.
    
    This function performs three main steps:
    1. Checks for existing genome files
    2. Downloads missing genomes using ncbi-genome-download
    3. Verifies and saves the downloaded genomes
    
    Args:
        genera: Set of bacterial genus names to download
        output_dir: Directory to save downloaded genomes
        config: Configuration parameters for download
        
    Returns:
        Dictionary mapping genus names to their genome file paths
        
    Raises:
        subprocess.CalledProcessError: If ncbi-genome-download fails
        FileNotFoundError: If no genomes are found in temp directory
    """

    output_dir.mkdir(exist_ok=True, parents=True)
    results = {}
    genera_to_download = set()

    # Stage 1: Check for existing genomes
    for genus in genera:
        if existing := check_existing_genome(genus, output_dir):
            if verify_genome_file(existing, genus):
                results[genus] = existing
                logger.info(f"Using existing genome for {genus}")
                continue 
        genera_to_download.add(genus)
    
    if not genera_to_download:
        return results 
    
    # Stage 2:  Downloads missing genomes using ncbi-genome-download
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            cmd = [
                    "ncbi-genome-download",
                    "bacteria",
                    "-F", config.file_format,
                    "--genera", ",".join(genera_to_download),
                    "--assembly-levels", config.assembly_level,
                    "--refseq-categories", config.refseq_category,
                    "--output-folder", temp_dir,
                    "--parallel", str(config.parallel),
                    "--retries", str(config.retries),
                    "--verbose"
                ]
            
            logger.info(f"Downloading genomes for: {', '.join(genera_to_download)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.debug(result.stdout)

            # Stage 3: Process downloaded files
            temp_genome_dir = Path(temp_dir) / "refseq" / "bacteria"
            if not temp_genome_dir.exists():
                raise FileNotFoundError(f"No genomes found in {temp_genome_dir}")
            
            for genus in genera_to_download:
                for genome_file in temp_genome_dir.glob("*/*_genomic.fna.gz"):
                    if verify_genome_file(genome_file, genus):
                        dest_file = output_dir / f"{genus}_reference.fna.gz"
                        dest_file.write_bytes(genome_file.read_bytes())
                        results[genus] = dest_file
                        logger.info(f"Successfully saved genome for {genus} to {dest_file}")
                        break
                else:
                    logger.warning(f"No valid genome found for {genus}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Download failed: {e.stderr}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")

    return results
        




        


