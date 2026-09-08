# CHTC high-dimensional experiments for the Subtrace/TempoSub paper

This directory contains the code and archived outputs for running the
100-biomarker synthetic baseline experiments on the Center for High Throughput
Computing (CHTC). The evaluated algorithms are BEBMS (with and without subtype
labels), SuStaIn GMM, and SuStaIn KDE.

## Exp5 and Exp8 correction

The original high-dimensional data were generated with the same 820-line
`sim_engine.py` used by the low-dimensional CHTC experiments. In that version,
low-dimensional biomarkers always used progression direction `+1`, but each
high-dimensional biomarker was assigned a fixed random `+1` or `-1` direction.
The direction flag is used only by the sigmoid measurement model, so this
affected only:

- Exp5: `sn_kjContinuousBeta_sigmoid`
- Exp8: `xiNearNormalWithNoise_kjContinuousBeta_sigmoid`

The intended convention was to let the biomarker parameters encode progression
direction and therefore use `+1` for every biomarker. Exp5 and Exp8 were
regenerated and rerun using this corrected convention. The other seven
experiments were unaffected and were not rerun.

During the redo, `config.yaml` enabled only Exp5 and Exp8. Consequently, the
CHTC staging input tarball and the regenerated `true_order_and_stages.json`
contained only the 100 replacement datasets (50 per experiment). This did not
delete the other experiments' result JSON files: the deletion commands in
`run.sh` were commented out, so the returned Exp5/8 files replaced matching old
files while the other results remained in `algo_results/`.

The current `algo_results/` therefore contains:

- corrected, rerun results for Exp5 and Exp8; and
- original results for Exp1–4, Exp6, Exp7, and Exp9.

All 400 current Exp5/8 algorithm JSON files match the replacement
`all_results.csv`. The 1,378 retained algorithm JSON files for the other
experiments match `all_results_before_redo_exp5_8.csv`.

## Result files

| File | Meaning |
| --- | --- |
| `all_results_before_redo_exp5_8.csv` | Historical results before the Exp5/8 correction. Its Exp5/8 rows are obsolete; its other rows are retained. |
| `all_results.csv` | Results produced after rerunning corrected Exp5/8. It contains only Exp5/8. |
| `final_all_results.csv` | Combined results created by `merge_two_results.py`; this is the file to use after regenerating it with the corrected merge script. |
| `algo_results/` | Raw result JSON files: corrected Exp5/8 plus retained results for the other experiments. |
| `true_order_and_stages.json` | Ground truth for the regenerated Exp5/8 datasets only. |

Because only two experiment names were enabled when `save_csv.py` created
`all_results.csv`, its positional numbering incorrectly assigned `E_Num = 1`
and `E_Num = 2` to Exp5 and Exp8. `merge_two_results.py` now corrects these to
5 and 8 using the experiment-name column before merging. It also removes the
old rows by experiment name rather than relying on `E_Num`.

To create the corrected combined CSV, run this command from this directory:

```sh
python3 merge_two_results.py
```

The expected output contains 2,227 rows with this experiment mapping:

| E_Num | Experiment | Rows |
| ---: | --- | ---: |
| 1 | `sn_kjOrdinalDM_xnjNormal` | 250 |
| 2 | `sn_kjOrdinalDM_xnjNonNormal` | 235 |
| 3 | `sn_kjOrdinalUniform_xnjNormal` | 250 |
| 4 | `sn_kjOrdinalUniform_xnjNonNormal` | 248 |
| 5 | `sn_kjContinuousBeta_sigmoid` | 250 |
| 6 | `sn_kjContinuousBeta_xnjNormal` | 250 |
| 7 | `sn_kjContinuousBeta_xnjNonNormal` | 244 |
| 8 | `xiNearNormalWithNoise_kjContinuousBeta_sigmoid` | 250 |
| 9 | `xiNearNormalWithNoise_kjContinuousBeta_xnjNormal` | 250 |

The total is below the theoretical 2,250 rows because 22 older algorithm jobs
are missing in Exp2, Exp4, and Exp7; one missing BEBMS result also means that
its corresponding Random Guessing row could not be generated.

## Archived CHTC workflow

- `gen.py` generates the configured datasets and ground-truth metadata.
- `gen.sh` historically regenerated the data, replaced the CHTC staging
  tarball, and submitted the jobs. It is destructive to the local `data/` and
  ground-truth files and should not be run merely to rebuild result tables.
- `gen_combo.py` creates `all_combinations.txt` from `config.yaml`.
- `run.sh` prepares and submits the HTCondor jobs.
- `run_mlhc.sub` defines file transfer and resource requests.
- `run_mlhc.sh` extracts the environment and staged data on each worker.
- `run_mlhc.py` runs all four algorithms for one dataset.
- `save_csv.py` converts raw result JSON files into a result table when the
  corresponding ground-truth metadata are available.
- `merge_two_results.py` combines the retained original rows with corrected
  Exp5/8 rows.

The CHTC staging area referenced by these scripts is historical and is not part
of this local directory.
