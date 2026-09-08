"""Merge the retained results with the corrected Exp5 and Exp8 results.

``all_results_before_redo_exp5_8.csv`` contains the original results for all
nine experiments. ``all_results.csv`` contains the rerun results for Exp5 and
Exp8, generated after fixing the high-dimensional sigmoid direction at +1.

Because only Exp5 and Exp8 were enabled when ``all_results.csv`` was created,
``save_csv.py`` numbered them 1 and 2 based on their positions in the shortened
configuration. Correct those numbers here before replacing the old Exp5/8 rows.
"""

import pandas as pd

ORIGINAL_FILE = "all_results_before_redo_exp5_8.csv"
APPEND_FILE = "all_results.csv"
OUTPUT_FILE = "final_all_results.csv"

EXP_NUMBER_MAP = {
    "sn_kjContinuousBeta_sigmoid": 5,
    "xiNearNormalWithNoise_kjContinuousBeta_sigmoid": 8,
}

df_original = pd.read_csv(ORIGINAL_FILE)
df_append = pd.read_csv(APPEND_FILE)

# The replacement CSV should contain only Exp5 and Exp8.
unexpected_experiments = set(df_append["E"].dropna().unique()) - set(EXP_NUMBER_MAP)
if unexpected_experiments:
    raise ValueError(
        "all_results.csv contains unexpected experiments: "
        f"{sorted(unexpected_experiments)}"
    )

# Assign the true experiment numbers explicitly. Do not rely on the shortened
# EXPERIMENT_NAMES list that was active when save_csv.py generated this CSV.
df_append = df_append.copy()
df_append["E_Num"] = df_append["E"].map(EXP_NUMBER_MAP)
if df_append["E_Num"].isna().any():
    raise ValueError("Some replacement rows could not be assigned an E_Num")
df_append["E_Num"] = df_append["E_Num"].astype(int)

# Remove the old Exp5/8 rows by experiment name, then append their rerun rows.
# Filtering by E is safer than trusting a potentially incorrect E_Num column.
df_filtered = df_original[~df_original["E"].isin(EXP_NUMBER_MAP)]
df_final = pd.concat([df_filtered, df_append], ignore_index=True)

df_final.to_csv(OUTPUT_FILE, index=False)

print(f"Done. Final row count: {len(df_final)}")
print(df_final.groupby(["E_Num", "E"]).size().to_string())
