#!/bin/bash

# directory
meta_dir="meta_dsk/txt_results"
dir_out="meta_genome_dsk"
genome_kmers="data/data_ggallus/genome/dsk/gal_genome.txt"
min_count=2  # min freq k-mer

# create output dir
mkdir -p "$dir_out"

# function for one file
process_file() {
    local meta_file="$1"
    local sample_name=$(basename "$meta_file" ".txt")
    local output_file="$dir_out/${sample_name}_shared_kmers.txt"
    local tmp_genome="$dir_out/genome_kmers_clean.tmp"
    local tmp_meta="$dir_out/meta_kmers_clean.tmp"

    # 1. Preparing temporary files (only if genome is not processed)
    [ -f "$tmp_genome" ] || awk '{print $1}' "$genome_kmers" > "$tmp_genome"

    # 2. Metagenome treatment
    awk '{print $1}' "$meta_file" > "$tmp_meta"

    # 3. Search for intersections
    if grep -Fxf "$tmp_genome" "$tmp_meta" > "$output_file"; then
        # 4. Frequency filtering
        if [ "$min_count" -gt 1 ]; then
            awk -v min="$min_count" '
                NR==FNR {genome[$1]=$2; next}
                $1 in genome && genome[$1] >= min {print $1, genome[$1]}
            ' "$genome_kmers" <(grep -Fwf "$output_file" "$meta_file") > "${output_file}.tmp" && \
            mv "${output_file}.tmp" "$output_file"
        fi

        # Statistic
        local count=$(wc -l < "$output_file")
        echo "${sample_name}: $count general k-mer (freq ≥$min_count)" >> "${dir_out}/summary.txt"
        return 0
    else
        echo "Processing error $meta_file" >> "${dir_out}/errors.log"
        return 1
    fi
}

# Exporting a function for GNU Parallel
export -f process_file
export genome_kmers dir_out min_count

# Basic processing cycle
for meta_file in "$meta_dir"/*_k15.h5.txt; do
    # retried file processing
    for attempt in {1..3}; do
        if process_file "$meta_file"; then
            break  # move on to the next file
        else
            sleep 10  #  before retrying
        fi
    done
done

# clean (only if successful)
[ $? -eq 0 ] && rm "$dir_out"/*.tmp 2>/dev/null
