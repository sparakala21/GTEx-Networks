import pandas as pd
import sqlite3

module_age_adjusted = pd.read_csv("GTEx_module_age_adjusted.tsv", sep = '\t')
module_bigtable_age_adjusted = pd.read_csv("GTEx_module_bigtable_age_adjusted.tsv", sep = '\t')
module_REACTOME_age_adjusted = pd.read_csv("GTEx_module_REACTOME_age_adjusted.tsv", sep = '\t')
network_age_adjusted = pd.read_csv("GTEx_network_age_adjusted.tsv", sep = '\t')
conn = sqlite3.connect('GTEx.db')

# Load DataFrame to SQLite table
module_age_adjusted.to_sql('GTEx_module_age_adjusted', conn, if_exists='replace', index=False)
module_bigtable_age_adjusted.to_sql('GTEx_module_bigtable_age_adjusted', conn, if_exists='replace', index=False)
module_REACTOME_age_adjusted.to_sql('GTEx_module_REACTOME_age_adjusted', conn, if_exists='replace', index=False)
network_age_adjusted.to_sql("GTEx_network_age_adjusted", conn, if_exists='replace', index=False)

# Close connection
conn.close()