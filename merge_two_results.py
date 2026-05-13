'''
Merge all_results_before_redo_exp5_8.csv and all_results. 

So right now, the data in the staging contian only exp5&8. This is because I changed the `sim_engine.py` 
which now uses fixed +-1 for sigmoid experiments. 

all_results_before_redo_exp5_8.csv contain all the results before I do this. 

all_results.csv now is only exp5&8 results. 

Note that in algo_results, I still have all the results. So the results in algo_results are complete.

I will rerun later, but not now. 

What I will do here is that i will delete all rows in all_results_before_redo_exp5_8.csv
whose E_num is 5. and just append all_results.csv to the end of the modified all_results_before_redo_exp5_8.csv.
'''

import pandas as pd

# File paths
original_file = "all_results_before_redo_exp5_8.csv"
append_file = "all_results.csv"

# Load CSVs
df_original = pd.read_csv(original_file)
df_append = pd.read_csv(append_file)

# Remove rows where E_Num is 5 or 8
df_filtered = df_original[~df_original["E_Num"].isin([5, 8])]

# Append new results
df_final = pd.concat([df_filtered, df_append], ignore_index=True)

# Save back to the original file
df_final.to_csv('final_all_results.csv', index=False)

print(f"Done. Final row count: {len(df_final)}")
