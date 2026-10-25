"""
05_sensitivity_and_baselines.py
===============================
Sensitivity checks and baselines reported in Sects. 4.7 and 5.2:

  * a conventional station-level QC baseline run on the same synthetic trials
  * POD counting only newly raised flags (vs flags already in the clean record)
  * sensitivity of attribution to the T4 weights, and share-vs-accuracy check
  * sensitivity of T6 / N5 to the 0.5 % tolerance, 0.70 threshold and window

The trial loop replays the random sequence of 03_error_injection.py, so the
4327 trials are the same ones.

Usage:  python src/05_sensitivity_and_baselines.py   (about 10 min)
Output: data/processed/sensitivity_trials.csv, data/processed/sensitivity_summary.json
"""
import numpy as np, pandas as pd, ncqc, json
from paths import PROCESSED

df = pd.read_csv(PROCESSED / 'merged.csv', parse_dates=['date']).set_index('date')
d0 = ncqc.build(df)
CLEAN = (d0.yr <= 2024) & (~d0.yr.isin([2016, 2018, 2019]))
par = ncqc.fit_tolerance(d0[CLEAN]); REF = CLEAN
ST = ['FS', 'SP', 'LR']
WEIGHTS = [0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]

# ---------- station-level baselines (fitted once on the clean record) ----------
def clim_params(d):
    p = {}
    for s in ST:
        lq = np.log(d0['Q_' + s])
        g = lq[CLEAN].groupby(d0.m[CLEAN])
        p[s] = (g.median(), g.apply(ncqc.mad))
    return p
CP = clim_params(d0)

def station_scores(d):
    """|z| of climatological test and |T| of temporal-continuity test, per station."""
    zc, zt = {}, {}
    T = ncqc.temporal_score(d, REF)
    for s in ST:
        lq = np.log(d['Q_' + s]); med, sc = CP[s]
        zc[s] = ((lq - d.m.map(med)) / d.m.map(sc)).abs()
        zt[s] = T['T_' + s].abs()
    flat = pd.concat([(d['Q_' + s].diff() == 0) for s in ST], axis=1).any(axis=1)
    return pd.DataFrame(zc).max(axis=1), pd.DataFrame(zt).max(axis=1).fillna(0), flat

def station_flag(d, thr):
    zc, zt, flat = station_scores(d)
    return (zc > thr) | (zt > thr) | flat

# calibrate the common threshold so that the station-level background flag
# rate on the clean record equals that of NCQC at |Z|>2
def detect(d, level=2.0):
    V = ncqc.T1(d); Z, _ = ncqc.closure_Z(d, par)
    A = ncqc.anomalies(d, REF); T = ncqc.temporal_score(d, REF); S = ncqc.spatial_score(d, A, REF)
    flag = (Z.abs() > level) | V['any']
    return flag, Z, T, S, V

base_flag, baseZ, _, _, _ = detect(d0)
FAR_NCQC = float(base_flag[CLEAN].mean())
zc0, zt0, fl0 = station_scores(d0)
def far_station(thr): return float(((zc0 > thr) | (zt0 > thr) | fl0)[CLEAN].mean())
grid = np.arange(1.5, 6.0, 0.01)
THR = float(grid[np.argmin([abs(far_station(t) - FAR_NCQC) for t in grid])])
# the individual tests alone at the same background rate
THR_C = float(grid[np.argmin([abs(float((zc0 > t)[CLEAN].mean()) - FAR_NCQC) for t in grid])])
THR_T = float(grid[np.argmin([abs(float((zt0 > t)[CLEAN].mean()) - FAR_NCQC) for t in grid])])
base_st = station_flag(d0, THR)
print('NCQC FAR', FAR_NCQC, 'station THR', THR, 'FAR', far_station(THR), 'THR_C', THR_C, 'THR_T', THR_T)

def attribute_w(T, S, V, w):
    sc = {s: w * T['T_' + s].abs() + (1 - w) * S['S_' + s].abs() for s in ST}
    SC = pd.DataFrame(sc)
    SC.loc[V['order_SP'], 'SP'] += 10; SC.loc[V['order_LR'], 'LR'] += 10
    return SC.fillna(0.0)

rng = np.random.default_rng(7)
clean_idx = np.where(CLEAN.values)[0]
rows = []
for typ in ['bias', 'spike', 'step', 'frozen']:
    for delta in [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]:
        for trial in range(160):
            s = rng.choice(ST)
            L = {'bias': rng.integers(2, 7), 'spike': 1, 'step': rng.integers(6, 13), 'frozen': rng.integers(1, 4)}[typ]
            i0 = rng.choice(clean_idx[(clean_idx > 6) & (clean_idx < len(d0) - int(L) - 6)])
            sl = slice(i0, i0 + int(L))
            if not CLEAN.values[sl].all(): continue
            sign = rng.choice([-1, 1])
            d = d0.copy(); col = 'Q_' + s
            if typ == 'frozen': d.iloc[sl, d.columns.get_loc(col)] = d[col].iloc[i0 - 1]
            else: d.iloc[sl, d.columns.get_loc(col)] = d[col].iloc[sl] * (1 + sign * delta)
            d = ncqc.build(d.drop(columns=['Qsum', 'res', 'rel', 'ndvi_bar', 'dndvi', 'regime', 'm', 'yr']))
            f, Z, T, S, V = detect(d)
            w_sl = f.iloc[sl]
            hit = bool(w_sl.any())
            pre = bool(base_flag.iloc[sl].any())
            new = bool((f & ~base_flag).iloc[sl].any())
            fs = station_flag(d, THR)
            zc, zt, flat = station_scores(d)
            r = dict(typ=typ, delta=delta, station=s, L=int(L), regime=d0.regime.iloc[sl].mode()[0],
                     hit=hit, pre=pre, new=new,
                     st_hit=bool(fs.iloc[sl].any()), st_new=bool((fs & ~base_st).iloc[sl].any()),
                     clim_hit=bool((zc > THR_C).iloc[sl].any()), temp_hit=bool((zt > THR_T).iloc[sl].any()),
                     flat_hit=bool(flat.iloc[sl].any()))
            if hit:
                for w in WEIGHTS:
                    SC = attribute_w(T, S, V, w); tot = SC.iloc[sl].sum()
                    r[f'att_{w}'] = tot.idxmax() == s
                    if w == 0.6: r['conf'] = float((tot / tot.sum()).max())
            rows.append(r)
R = pd.DataFrame(rows); R.to_csv(PROCESSED / 'sensitivity_trials.csv', index=False)
print('trials', len(R))

out = {'FAR_NCQC': FAR_NCQC, 'THR_station': THR, 'THR_clim': THR_C, 'THR_temp': THR_T,
       'FAR_station_combined': far_station(THR)}
pod = R.groupby(['typ', 'delta'])[['hit', 'new', 'pre', 'st_hit', 'st_new', 'clim_hit', 'temp_hit', 'flat_hit']].mean().round(3)
print(pod.to_string())
out['pod'] = {f'{a}|{b}': v for (a, b), v in pod.to_dict('index').items()}
H = R[R.hit]
att = {}
for w in WEIGHTS:
    att[w] = H.groupby('station')[f'att_{w}'].mean().round(3).to_dict()
    att[w]['all'] = round(float(H[f'att_{w}'].mean()), 3)
print(pd.DataFrame(att).T.to_string())
out['att_weights'] = {str(k): v for k, v in att.items()}
H = H.assign(cbin=pd.cut(H.conf, [0, 0.4, 0.5, 0.6, 1.0]))
cal = H.groupby('cbin', observed=True)['att_0.6'].agg(['mean', 'size']).round(3)
print(cal.to_string()); out['calibration'] = {str(k): v for k, v in cal.to_dict('index').items()}
cal_st = H.groupby(['station', 'cbin'], observed=True)['att_0.6'].agg(['mean', 'size']).round(3)
print(cal_st.to_string())

# ---------- episode-level weight sensitivity on the real record ----------
res = ncqc.run(d0, REF)
ep_rows = []
for w in WEIGHTS:
    SC = attribute_w(res['T'], res['S'], res['V'], w)
    E = ncqc.attribute_episodes(d0, SC, res['flag'], res['Z'])
    E = E[E.n > 1]
    for _, e in E.iterrows():
        ep_rows.append(dict(w=w, start=str(e.start.date()), source=e.source, conf=e.confidence))
EP = pd.DataFrame(ep_rows).pivot(index='start', columns='w', values='source')
EPc = pd.DataFrame(ep_rows).pivot(index='start', columns='w', values='conf')
print(EP.to_string()); print(EPc.to_string())
out['episode_sources'] = EP.to_dict(); out['episode_conf'] = EPc.to_dict()

# ---------- T6 sensitivity ----------
t6 = []
ice = d0[d0.regime == 'Ice-affected']
for tight in [0.25, 0.5, 0.75, 1.0]:
    for win in [8, 12, 16]:
        frac = ncqc.redundancy(d0, win=win, tight=tight)
        for thr in [0.6, 0.7, 0.8]:
            above = frac[frac > thr]
            yrs = sorted(set(above.index.year))
            in_2020 = frac[(frac.index.year >= 2020) & (frac.index.year <= 2024)]
            pre = frac[frac.index.year <= 2019]
            t6.append(dict(tight=tight, win=win, thr=thr,
                           fires_2020_24=bool((in_2020 > thr).any()),
                           fires_before_2020=bool((pre > thr).any()),
                           fires_2025_26=bool((frac[frac.index.year >= 2025] > thr).any()),
                           max_pre2020=round(float(pre.max()), 2), max_2020_24=round(float(in_2020.max()), 2),
                           years=','.join(map(str, yrs))))
T6 = pd.DataFrame(t6); print(T6.to_string()); out['t6'] = T6.to_dict('records')
# fraction of tight months by epoch for several tolerances
ep = pd.cut(ice.yr, [1999, 2019, 2024, 2026], labels=['2000-19', '2020-24', '2025-26'])
out['tight_by_epoch'] = {str(t): (ice.rel.abs() < t).groupby(ep, observed=True).mean().round(2).to_dict() for t in [0.25, 0.5, 0.75, 1.0]}
print(out['tight_by_epoch'])
json.dump(out, open(PROCESSED / 'sensitivity_summary.json', 'w'), indent=1, default=str)
