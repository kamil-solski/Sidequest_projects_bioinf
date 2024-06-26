"""
 Poniższy program parsuje tabele uzyskaną za pomocą mmseqs2 (w formacie tsv) oraz mapuje sekwencje do GO. Na koniec zwraca plik z wynikami w csv
"""

import pandas as pd
import os
import sys

# Parsowanie pliku pfam2go 
def parse_pfam2go(pfam2go_file):
    pfam2go = {}
    with open(pfam2go_file, "r") as f:
        for line in f.readlines()[6:]:
            try:
                pfam, prot = line.split(" > ")[0].split()[0].split(":")[1].strip(), line.split(" > ")[0].split()[1].strip()  # wyciąga pierwszą część linijki np. Pfam:PF0001 7tm_1 (line.split(" > ")[0]), rozszczepia na ['Pfam:PF00001', '7tm_1'] (line.split(" > ")[0].split()), na koniec pobiera PF00001 jako pfam i 7tm_1 jako prot
                fun, goid = line.split(" > ")[1].split(";")[0].strip(), line.split(" > ")[1].split(";")[1].strip()  # wyciąga drugą część linijki np. GO:G protein-coupled receptor activity ; GO:0004930, rozszczepia na ['GO:G protein-coupled receptor activity', 'GO:0004930'] (line.split(" > ")[1].split(";")), na koniec pobiera 'G protein-coupled receptor activity' jako fun i 'GO:0004930' jako goid
            except Exception as e:  # pomija potencjalne błędy w linijkach pliku pfam2go.txt
                print(f"Error parsing line: {line}\n{e}")
                continue
            pfam2go.setdefault(pfam, []).append((pfam, goid, prot, fun))
    return pfam2go  # otrzymamy słownik: {pfam_id: [go_terms]}
  

# Mapowanie sekwencji do GO
def map_sequences_to_go(df, pfam2go):
    not_found = 0  # inicjalizacja zmiennej zliczającej ile razy dany Pfam ID w pliku z sekwencjami nie został znaleziony w pfam2go
    df["pfam"] = [x.split(".")[0] for x in df["PFAM_model"]]  # tworzymy nową kolumne z samym identyfikatorem Pfam
    in_pfam = df["pfam"]  # tworzenie listy aby móc po niej iterować
    sequence_to_go = []
    for pf in in_pfam:
        if pf in pfam2go:  # sprawdzanie czy dany Pfam ID istnieje w słowniku pfam2go. Jeżeli istnieje to iteruje po liście GO terms (przykład: GO:0004930, GO:0007186...)
            for g in pfam2go[pf]:
                sequence_to_go.append(g)  # na koniec dodaje GO terms do listy
        else:
            not_found += 1  # zlicza przypadki w których dany Pfam ID nie został znaleziony w pfam2go
    return sequence_to_go, not_found
  
  
def calculate_go(sequence_to_go):  # zlicza powtórzenia dla konkretnych GO terms w liście sequence_to_go
    counts = {}  # tylko unikalne
    total = 0  # wszytkie
    for g in sequence_to_go:
        if g[1] not in counts:
            counts[g[1]] = 0
        counts[g[1]] += 1
        total += 1
    return counts, total

def save_to_csv(counts, total, output_file):
    data = [(goid, count/total) for goid, count in counts.items()]
    df = pd.DataFrame(data, columns=["GO ID", "Frequency"])
    df.to_csv(output_file, index=False)
  

# Wykonanie wszytkich funkcji
def main():
    if len(sys.argv) != 4:
        print("Wrong arguments. Pfam.py <input_tsv_file> <pfam2go_file> <output_csv_file>")
        sys.exit(1)

    input_file = sys.argv[1]  # podaj plik z sekwencjami
    pfam2go_file = sys.argv[2]  # podaj plik pfam2go
    output_file = sys.argv[3]  # podaj nazwę pliku z wynikami w csv

    # Wczytanie pliku tsv
    df = pd.read_csv(input_file, sep='\t', header=None, 
                     names=['PFAM_model', 'Sequence_ID', 'alnScore', 'seqIdentity', 'eVal', 
                            'qStart', 'qEnd', 'qLen', 'tStart', 'tEnd', 'tLen'])

    pfam2go = parse_pfam2go(pfam2go_file)

    sequences_to_go, not_found = map_sequences_to_go(df, pfam2go)

    counts, total = calculate_go(sequences_to_go)

    # Wypisz w terminalu wszystkie wartości GO
    for goid, count in counts.items():
        print(f"{goid},{count/total:.6f}")
    
    # Zapisz w pliku csv
    save_to_csv(counts, total, output_file)
        
if __name__ == "__main__":
    main()