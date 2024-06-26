#!/bin/bash

# Sprawdzanie plików
files=("Pfam.py" "pfam2go.txt" "FunctionalAnnotationContigs.zip")

all_present=true

for file in "${files[@]}"; do
  if [ ! -e "$file" ]; then
    echo "Error: $file is not present in folder."
    all_present=false
    
  fi
done

if [ "$all_present" = false ]; then
  echo "Stopping program due to missing file/s."
  exit 1
fi

echo "All files present!"

# tworzenie środowiska
if ! conda env list | grep -q "\<Metagenomic_3\>"; then
    echo "Creating conda environment: Metagenomic_3"
    conda create -n Metagenomic_3 python=3.8
fi

conda activate Metagenomic_3

# isntalowanie potrzebnych paczek
conda install bioconda::mmseqs2
conda install bioconda::prodigal
pip install plotly
conda install conda-forge::dash-bio

unzip FunctionalAnnotationContigs.zip -d ./FunctionalAnnotationContigs

# wykonanie prodigal
cd FunctionalAnnotationContigs

prodigal -i Finland_DrinkingWater.ERZ1758444_FASTA.fasta -o finland_prodigal -a finland_trans.faa -p meta & pid=$!

prodigal -i Germany_WasteWater.ERZ3455850_FASTA.fasta -o germany_prodigal -a germany_trans.faa -p meta

wait $pid

# klastrowanie i przeszukiwanie za pomocą mmsqs2
mmseqs databases Pfam-A.seed Pfam-A_seed tmp

mmseqs createdb finland_trans.faa finland.mmseqs
mmseqs search --max-seqs 30000000 Pfam-A_seed finland.mmseqs pfam_finland_bin_out tmp/ -s 7 --threads 6
mmseqs createtsv Pfam-A_seed finland.mmseqs pfam_finland_bin_out finland.tsv

mmseqs createdb germany_trans.faa germany.mmseqs
mmseqs search --max-seqs 30000000 Pfam-A_seed germany.mmseqs pfam_germany_bin_out tmp/ -s 7 --threads 6
mmseqs createtsv Pfam-A_seed germany.mmseqs pfam_germany_bin_out germany.tsv

# Parsowanie otrzymanej tabeli tsv (przypisanie odpowiednich nazw kolumn) oraz zmapowanie pfam2go.txt z sekwencją. Skrypt przesyłam w osobnym pliku Pfam.py
cd ..

# podaję nazwę skryptu nazwę pliku z sekwencjami, pfam2go i nazwy outputu w formacie csv 
python Pfam.py FunctionalAnnotationContigs/finland.tsv pfam2go.txt fin.csv

python Pfam.py FunctionalAnnotationContigs/germany.tsv pfam2go.txt ger.csv


# wizualizacja clustergram obu outputów
python Visualize.py fin.csv ger.csv
