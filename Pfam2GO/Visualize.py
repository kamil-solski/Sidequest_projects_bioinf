import sys
import pandas as pd
import plotly
import dash_bio

try:
    fin = pd.read_csv(sys.argv[1], skiprows=1, names=["go", "count_fin"])  # skiprows=1 ponieważ plik zawiera nazwę kolumny w pierwszym wierszu
    ger = pd.read_csv(sys.argv[2], skiprows=1, names=["go", "count_ger"])
# zatrzymuje wykonywanie jeżeli został wczytany zły plik
except Exception as e:  
    print(f"Error reading files: {e}")
    sys.exit(1)

if not pd.api.types.is_numeric_dtype(fin['count_fin']) or not pd.api.types.is_numeric_dtype(ger['count_ger']):
    print("Count columns must be numeric.")
    sys.exit(1)  # w razie błędnych danych (nie numerycznych) w zliczeniach program zatrzymuje wykonywanie. Jest to mało prawdopodobne 
 

combined = fin.merge(ger, how="outer", on="go")  # połączenie danych z finlandi i niemiec do jednej tabeli
combined["count_ger"] = combined["count_ger"].fillna(0)  # zastępuje puste wartości zerami
combined["count_fin"] = combined["count_fin"].fillna(0)
combined["mean"] = (combined["count_ger"] + combined["count_fin"]) / 2  # dodaje kolumnę średnia
combined.sort_values(by="mean", ascending=False, inplace=True)

# wyciągniecie tylko 100 najwyzszych wyników względem średniej. Jeżeli danych jest mniej niż 100 wtedy weźmie tą mniejszą liczbę
top_n = min(len(combined), 100)  
to_vis = combined.head(top_n)
to_vis2 = to_vis[["count_ger", "count_fin"]]
ll = list(to_vis["go"])

fig = dash_bio.Clustergram(data=to_vis2, column_labels=["count_ger", "count_fin"], row_labels=ll, height=1600, width=700)
fig.write_html("Clustergram.html")
