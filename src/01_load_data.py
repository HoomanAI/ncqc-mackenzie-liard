"""
01_load_data.py
===============
Reads the three station workbooks (WSC discharge plus ERA5-Land and MODIS
values extracted with Google Earth Engine) for the Mackenzie-Liard confluence
and writes one tidy monthly table, data/processed/merged.csv, which every
later script reads.

Stations
    FS  10GC001  Mackenzie River at Fort Simpson   (closure node, downstream)
    SP  10FB006  Mackenzie River at Strong Point   (upstream main-stem node)
    LR  10ED002  Liard River near the mouth        (upstream tributary node)

Each workbook holds, at monthly resolution for 2000-01 to 2026-07:
    Runoff (m)          ERA5-Land surface runoff at the gauge pixel
    Precipitation (m)   ERA5-Land total precipitation at the gauge pixel
    Q (m3/s)            reported monthly mean discharge
    NDVI                MODIS 16-day NDVI composited to monthly, gauge pixel

Usage:  python src/01_load_data.py
Output: data/processed/merged.csv
"""

import pandas as pd, numpy as np
from paths import RAW, PROCESSED

FILES = {
 'FS': (RAW / 'Fort_Simpson_GEE_2000-2026.xlsx', '10GC001'),
 'LR': (RAW / 'Liard_River_GEE_2000-2026.xlsx', '10ED002'),
 'SP': (RAW / 'Strong_Point_GEE_2000-2026.xlsx', '10FB006'),
}

def load(p):
    raw = pd.read_excel(p, header=None)
    df = raw.iloc[2:,0:5].copy()
    df.columns=['date','runoff','precip','Q','ndvi']
    df['date']=pd.to_datetime(df['date'])
    for c in ['runoff','precip','Q','ndvi']:
        df[c]=pd.to_numeric(df[c],errors='coerce')
    return df.set_index('date')

d={}
for k,(p,sid) in FILES.items():
    d[k]=load(p)

# the three workbooks must cover exactly the same months: assigning columns
# aligns on the index, so a mismatch would silently inject NaN instead of failing
ref=d['FS'].index
for k,v in d.items():
    if not v.index.equals(ref):
        raise ValueError(f'{k}: date index differs from FS '
                         f'({len(v.index)} vs {len(ref)} rows; '
                         f'first mismatch {v.index.symmetric_difference(ref)[:3].tolist()})')

df = pd.DataFrame({'date':ref})
df=df.set_index('date')
for k in d:
    for c in ['Q','runoff','precip','ndvi']:
        df[f'{c}_{k}']=d[k][c]
assert df.notna().all().all(), 'merged table contains NaN - check the workbooks'
df.to_csv(PROCESSED / 'merged.csv')
print(df.shape)
print(df.isna().sum())
print(df.describe().T.to_string())
