from scripts.download_metagenome import download_genomes
from pathlib import Path
import pandas as pd
import os


current_dir = os.path.dirname(os.path.dirname(os.getcwd()))
file_path = os.path.join(current_dir, 'species.csv')

genera_df = pd.read_csv(
    file_path,
    index_col=0, 
    skipinitialspace=True, 
    dtype={'Genus': str, 'CCA1': float, 'CCA2': float}
)


output_dir = Path("output_downloaded_metagenome")

genera_df['Genus'] = genera_df['Genus'].str.strip()

genera = genera_df['Genus'].tolist()

results = download_genomes(genera, output_dir)

print("Downloaded genomes:")
for genus, path in results.items():
    print(f"{genus}: {path}")




    
