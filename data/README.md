# Data

## raw/

One workbook per station, monthly values from January 2000 to July 2026 (319 months, no gaps).

| File | Station | WSC ID | Lat, Lon |
|---|---|---|---|
| `Fort_Simpson_GEE_2000-2026.xlsx` | Mackenzie River at Fort Simpson (downstream) | 10GC001 | 61.8685, −121.3589 |
| `Strong_Point_GEE_2000-2026.xlsx` | Mackenzie River at Strong Point (upstream main stem) | 10FB006 | 61.8165, −120.7917 |
| `Liard_River_GEE_2000-2026.xlsx` | Liard River near the mouth (upstream tributary) | 10ED002 | 61.7426, −121.2278 |

The first sheet of each workbook has a two-line header followed by these columns:

| Column | Units | Source |
|---|---|---|
| Date | first day of month | – |
| Runoff | m | ERA5-Land surface runoff at the grid cell containing the gauge |
| Precipitation | m | ERA5-Land total precipitation at the same cell |
| Q | m³ s⁻¹ | WSC monthly mean discharge |
| NDVI | – | MODIS MOD13Q1 NDVI at the gauge pixel, averaged to the month |

The `Qmax_daily` sheet lists daily peak flows for selected events. The analysis does not use it.

Discharge record status:

- **Jan 2000 – Dec 2024:** approved records (HYDAT, release of 17 July 2026), reported to about three significant digits.
- **Jan 2025 – Jul 2026:** provisional records from Wateroffice. Monthly values are means of the provisional daily values that were available. A few months were missing some days, and several winter values look like constant ice-period estimates. These records may still be revised by WSC.

## processed/

Written by the scripts in `src/`. They are committed so that results can be checked without re-running anything.

| File | Written by | Contents |
|---|---|---|
| `merged.csv` | 01 | All three stations side by side (`Q_`, `runoff_`, `precip_`, `ndvi_` + FS/SP/LR) |
| `injection.csv` | 03 | One row per synthetic trial: error type, size, station, hit, attribution |
| `episodes.csv` | 04 | Flagged episodes with attributed station, share and max \|Z\| |
| `summary.json` | 04 | Tolerance parameters, flag counts, POD tables, attribution rates |
| `sensitivity_trials.csv` | 05 | The same trials with baseline results and alternative attribution weights |
| `sensitivity_summary.json` | 05 | Baseline POD, weight sensitivity, T6 threshold sensitivity |

## Sources and terms

- Environment and Climate Change Canada, Water Survey of Canada: National Water Data Archive (HYDAT) and Wateroffice. Open Government Licence – Canada.
- Muñoz-Sabater, J. et al. (2021) ERA5-Land: a state-of-the-art global reanalysis dataset for land applications. *Earth Syst. Sci. Data*, 13, 4349–4383. Copernicus licence.
- Didan, K. (2015) MOD13Q1 MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V006. NASA LP DAAC, https://doi.org/10.5067/MODIS/MOD13Q1.006
