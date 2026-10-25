"""
03_error_injection.py
=====================
Skill assessment of the NCQC detector by synthetic error injection.

Corruptions applied to one randomly chosen station over a randomly chosen
window of the "clean" reference record (2000-2024, excluding the 2016 and
2018-2019 anomaly epochs):

    bias    multiplicative offset of +/- delta held for 2-6 months
    spike   single-month multiplicative offset of +/- delta
    step    multiplicative offset of +/- delta held for 6-12 months
    frozen  the previous month's value repeated for 1-3 months
            (delta is NOT applied to this type; it is kept only as a trial label
             so that the four types share one experimental grid, which is why
             Fig. 9a plots the frozen curve against a nominal abscissa)

Tolerance parameters are frozen on the clean record and are NOT refitted for
each trial, so the experiment measures out-of-sample detection skill.
4327 trials in total (160 per error type and magnitude).

Usage:  python src/03_error_injection.py   (needs data/processed/merged.csv)
Output: data/processed/injection.csv
"""

import numpy as np, pandas as pd, ncqc
from paths import PROCESSED

df=pd.read_csv(PROCESSED / 'merged.csv',parse_dates=['date']).set_index('date')
d0=ncqc.build(df)
CLEAN=(d0.yr<=2024)&(~d0.yr.isin([2016,2018,2019]))
par=ncqc.fit_tolerance(d0[CLEAN])          # frozen tolerance parameters
REF=CLEAN

def detect(d, level=2.0):
    V=ncqc.T1(d)
    Z,_=ncqc.closure_Z(d,par)
    A=ncqc.anomalies(d,REF); T=ncqc.temporal_score(d,REF); S=ncqc.spatial_score(d,A,REF)
    src,SC=ncqc.attribute(d,T,S,V)
    flag=(Z.abs()>level)|V['any']
    return flag,Z,SC,src

base_flag,baseZ,_,_=detect(d0)
# false-alarm rate on clean months
FAR0={}
for lev in [1.5,2,2.5,3,4]:
    f,_,_,_=detect(d0,lev)
    FAR0[lev]=float(f[CLEAN].mean())
print('Baseline flag rate on clean record:',{k:round(v,3) for k,v in FAR0.items()})

rng=np.random.default_rng(7)
clean_idx=np.where(CLEAN.values)[0]
types=['bias','spike','step','frozen']
rows=[]
LEVEL=2.0
for typ in types:
    for delta in [0.02,0.05,0.10,0.15,0.20,0.30,0.50]:
        for trial in range(160):
            s=rng.choice(['FS','SP','LR'])
            L={'bias':rng.integers(2,7),'spike':1,'step':rng.integers(6,13),'frozen':rng.integers(1,4)}[typ]
            i0=rng.choice(clean_idx[(clean_idx>6)&(clean_idx<len(d0)-int(L)-6)])
            sl=slice(i0,i0+int(L))
            if not CLEAN.values[sl].all(): continue
            sign=rng.choice([-1,1])
            d=d0.copy()
            col='Q_'+s
            if typ=='frozen':
                d.iloc[sl, d.columns.get_loc(col)]=d[col].iloc[i0-1]
            else:
                d.iloc[sl, d.columns.get_loc(col)]=d[col].iloc[sl]*(1+sign*delta)
            d=ncqc.build(d.drop(columns=['Qsum','res','rel','ndvi_bar','dndvi','regime','m','yr']))
            f,Z,SC,src=detect(d,LEVEL)
            hit=bool(f.iloc[sl].any())
            maxZ=float(Z.iloc[sl].abs().max()); viol=bool(ncqc.T1(d)['any'].iloc[sl].any())
            att=None
            if hit:
                tot=SC.iloc[sl].sum(); att=tot.idxmax(); conf=float((tot/tot.sum()).max())
            else: conf=np.nan
            rg=d0.regime.iloc[sl].mode()[0]
            # realised relative perturbation; for 'frozen' this is whatever the
            # repeated value implies, not the nominal delta
            eff=abs(d.iloc[sl][col].values/d0.iloc[sl][col].values-1).mean()
            rows.append(dict(typ=typ,delta=delta,station=s,L=int(L),regime=rg,hit=hit,
                             att=att,correct=(att==s),conf=conf,eff=eff,maxZ=maxZ,viol=viol))
R=pd.DataFrame(rows)
R.to_csv(PROCESSED / 'injection.csv',index=False)
print(R.groupby(['typ','delta'])['hit'].mean().round(2).unstack().to_string())
print()
print(R[R.typ.isin(['bias','spike','step'])].groupby(['regime','delta'])['hit'].mean().round(2).unstack().to_string())
print()
print(R[R.hit].groupby(['typ'])['correct'].mean().round(2).to_string())
print(R[R.hit].groupby(['typ','station'])['correct'].mean().round(2).unstack().to_string())
print()
print(R[R.typ=='frozen'].groupby('regime')['hit'].mean().round(2).to_string())
