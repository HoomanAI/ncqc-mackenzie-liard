# Changelog

## 1.1.0 – 2026-10-25

- Added `05_sensitivity_and_baselines.py`: station-level QC baseline, POD for newly raised flags only, attribution-weight sensitivity and T6 threshold sensitivity (Sects. 4.7 and 5.2 of the paper).
- Figure 3a now shades ice-affected months; Figure 6 was redrawn with sharper lines and larger type; the Figure 9d label for Strong Point now reads "main stem".
- Scripts now run from the repository root with a fixed `data/` and `figures/` layout (`src/paths.py`).
- No change to any computed result: `merged.csv`, `injection.csv`, `episodes.csv` and `summary.json` are byte-identical to 1.0.1.

## 1.0.1 – 2026-10-20

- T6 (`redundancy()`) now windows over consecutive ice-affected months only; previously it mixed in open-water months, so N5 could never trigger.
- `attribute_episodes()` returns the max |Z| per episode (Table 4).
- `01_load_data.py` checks that the three workbooks share the same dates.
- Removed a duplicated branch in `03_error_injection.py` (no effect on results).

## 1.0.0 – 2026-09-15

- First version, submitted with the manuscript.
