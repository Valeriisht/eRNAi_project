from download_metagenome import download_genomes
from pathlib import Path
import pandas as pd
import os

genera_df = pd.read_csv(os.path.dirname(os.path.dirname(os.getcwd())) + '/filtered_genera.csv')


output_dir = Path("output_downloaded_metagenome")

genera = genera_df["Genus"].tolist()

results = download_genomes(genera, output_dir)

print("Downloaded genomes:")
for genus, path in results.items():
    print(f"{genus}: {path}")


    
