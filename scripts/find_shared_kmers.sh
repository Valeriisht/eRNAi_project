#!/bin/bash

# Директории
meta_dir="meta_dsk/txt_results"
dir_out="meta_genome_dsk"
genome_kmers="data/data_ggallus/genome/dsk/gal_genome.txt"
min_count=2  # Минимальная частота k-mer

# Создаем выходную директорию
mkdir -p "$dir_out"

# Функция для обработки одного файла
process_file() {
    local meta_file="$1"
    local sample_name=$(basename "$meta_file" ".txt")
    local output_file="$dir_out/${sample_name}_shared_kmers.txt"
    local tmp_genome="$dir_out/genome_kmers_clean.tmp"
    local tmp_meta="$dir_out/meta_kmers_clean.tmp"

    # 1. Подготовка временных файлов (только если genome не обработан)
    [ -f "$tmp_genome" ] || awk '{print $1}' "$genome_kmers" > "$tmp_genome"

    # 2. Обработка метагенома
    awk '{print $1}' "$meta_file" > "$tmp_meta"

    # 3. Поиск пересечений
    if grep -Fxf "$tmp_genome" "$tmp_meta" > "$output_file"; then
        # 4. Фильтрация по частоте
        if [ "$min_count" -gt 1 ]; then
            awk -v min="$min_count" '
                NR==FNR {genome[$1]=$2; next}
                $1 in genome && genome[$1] >= min {print $1, genome[$1]}
            ' "$genome_kmers" <(grep -Fwf "$output_file" "$meta_file") > "${output_file}.tmp" && \
            mv "${output_file}.tmp" "$output_file"
        fi

        # Статистика
        local count=$(wc -l < "$output_file")
        echo "${sample_name}: $count общих k-mer (частота ≥$min_count)" >> "${dir_out}/summary.txt"
        return 0
    else
        echo "Ошибка обработки $meta_file" >> "${dir_out}/errors.log"
        return 1
    fi
}

# Экспорт функции для GNU Parallel
export -f process_file
export genome_kmers dir_out min_count

# Основной цикл обработки
for meta_file in "$meta_dir"/*_k15.h5.txt; do
    # Пытаемся обработать файл с повторными попытками
    for attempt in {1..3}; do
        if process_file "$meta_file"; then
            break  # Успешно - переходим к следующему файлу
        else
            sleep 10  # Ждем перед повторной попыткой
        fi
    done
done

# Очистка (только если все успешно)
[ $? -eq 0 ] && rm "$dir_out"/*.tmp 2>/dev/null