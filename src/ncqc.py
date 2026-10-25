"""
ncqc.py
=======
NCQC - Network-Constrained Quality Control for a hydrometric nodal control
volume.  Implements the seven-stage framework described in the manuscript:

    T0  provenance signature        sigdigits()
    T1  hard admissibility          T1()
    T2  regime classification       classify_regime()
    T3  nodal closure test          fit_tolerance(), closure_Z()
    T4  attribution                 temporal_score(), spatial_score(), attribute()
    T5  episode / change point      episodes(), attribute_episodes(), cusum()
    T6  independence diagnostic     redundancy()

run() chains T1-T4 and returns every intermediate diagnostic.
This module is imported by the numbered scripts; it has no side effects.
"""
import numpy as np, pandas as pd

MADK = 1.4826
def mad(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return MADK*np.median(np.abs(x-np.median(x))) if len(x) else np.nan

def classify_regime(ndvi_bar, dndvi):
    out=[]
    for n,d in zip(ndvi_bar,dndvi):
        if np.isfinite(d) and d>=0.15: out.append('Transition')
        elif n<0.15: out.append('Ice-affected')
        else: out.append('Open-water')
    return np.array(out)

def build(df):
    d=df.copy()
    d['Qsum']=d.Q_SP+d.Q_LR
    d['res']=d.Q_FS-d.Qsum
    d['rel']=d.res/d.Q_FS*100
    d['ndvi_bar']=d[['ndvi_FS','ndvi_LR','ndvi_SP']].mean(axis=1)
    d['dndvi']=d.ndvi_bar.diff().abs()
    d['regime']=classify_regime(d.ndvi_bar.values,d.dndvi.values)
    d['m']=d.index.month; d['yr']=d.index.year
    return d

# ---------------- T1 hard admissibility -------------------------------------
def T1(d, eps_neg=0.0):
    """Hard physical violations: non-positive flow, upstream>downstream,
    implied negative tributary contribution."""
    v = pd.DataFrame(index=d.index)
    v['nonpos'] = (d[['Q_FS','Q_SP','Q_LR']]<=0).any(axis=1)
    v['order_SP'] = d.Q_SP > d.Q_FS*(1+eps_neg)
    v['order_LR'] = d.Q_LR > d.Q_FS*(1+eps_neg)
    v['any'] = v[['nonpos','order_SP','order_LR']].any(axis=1)
    return v

# ---------------- T2/T3 regime-conditioned closure ---------------------------
def fit_tolerance(d_ref):
    """Robust, regime-specific proportional error scale sigma(Q)=k_r*Q + s0_r."""
    par={}
    for rg,g in d_ref.groupby('regime'):
        k = mad(g.res/g.Q_FS)          # proportional component
        s0 = mad(g.res)*0.25           # small additive floor
        b  = np.median(g.res/g.Q_FS)   # regime bias
        par[rg]=dict(k=max(k,1e-4), s0=s0, bias=b, n=len(g))
    return par

def closure_Z(d, par):
    k  = d.regime.map(lambda r: par[r]['k'])
    s0 = d.regime.map(lambda r: par[r]['s0'])
    b  = d.regime.map(lambda r: par[r]['bias'])
    sigma = k*d.Q_FS + s0
    return (d.res - b*d.Q_FS)/sigma, sigma

# ---------------- T4 temporal + spatial plausibility (attribution) -----------
def anomalies(d, ref_mask):
    a={}
    for s in ['FS','SP','LR']:
        clim = np.log(d['Q_'+s]).groupby(d.m).transform(lambda x: x[ref_mask.reindex(x.index,fill_value=False)].median()
                                                        if ref_mask.reindex(x.index,fill_value=False).any() else x.median())
        a['a_'+s]=np.log(d['Q_'+s])-clim
    return pd.DataFrame(a,index=d.index)

def temporal_score(d, ref_mask):
    """Local continuity: departure of log Q from the mean of its neighbours in time."""
    out={}
    for s in ['FS','SP','LR']:
        lq=np.log(d['Q_'+s])
        pred=0.5*(lq.shift(1)+lq.shift(-1))
        e=lq-pred
        sc=pd.Series(index=d.index,dtype=float)
        for rg in d.regime.unique():
            m=(d.regime==rg)
            s_rg=mad(e[m & ref_mask])
            sc[m]=e[m]/max(s_rg,1e-6)
        out['T_'+s]=sc
    return pd.DataFrame(out)

def spatial_score(d, A, ref_mask):
    """Departure from the inter-station anomaly regressions estimated on reference months."""
    out={}
    pairs={'FS':['SP','LR'],'SP':['FS','LR'],'LR':['FS','SP']}
    for s,others in pairs.items():
        X=A[['a_'+o for o in others]].values; y=A['a_'+s].values
        M=ref_mask.values & np.isfinite(X).all(axis=1) & np.isfinite(y)
        Xm=np.column_stack([np.ones(M.sum()),X[M]])
        beta,_,_,_=np.linalg.lstsq(Xm,y[M],rcond=None)
        pred=beta[0]+X@beta[1:]
        e=y-pred
        out['S_'+s]=e/max(mad(e[M]),1e-6)
    return pd.DataFrame(out,index=A.index)

def attribute(d, T, S, V):
    """Composite implausibility score per station; argmax = attributed source."""
    sc={}
    for s in ['FS','SP','LR']:
        base=0.6*T['T_'+s].abs()+0.4*S['S_'+s].abs()
        sc[s]=base
    SC=pd.DataFrame(sc)
    SC.loc[V['order_SP'],'SP']+=10
    SC.loc[V['order_LR'],'LR']+=10
    SCf=SC.fillna(0.0)
    return SCf.idxmax(axis=1), SC

# ---------------- T5 CUSUM change point --------------------------------------
def cusum(z, k=0.5):
    zp=np.nan_to_num(z.values); Sp=np.zeros(len(zp)); Sn=np.zeros(len(zp))
    for i in range(1,len(zp)):
        Sp[i]=max(0,Sp[i-1]+zp[i]-k)
        Sn[i]=min(0,Sn[i-1]+zp[i]+k)
    return pd.Series(Sp,index=z.index), pd.Series(Sn,index=z.index)

# ---------------- T6 redundancy / independence -------------------------------
def redundancy(d, win=12, tight=0.5, regime='Ice-affected'):
    """T6: fraction of ice-affected months, within a moving window of `win`
    consecutive ice-affected samples, whose relative closure error is
    implausibly small.  A sustained fraction above 0.7 raises flag N5: the
    agreement is too close for closure alone to count as independent evidence
    of accuracy (it does not prove that neighbour data were used).

    Note the window counts ICE-AFFECTED SAMPLES, not calendar months: with
    four to six ice months a year, win=12 spans about three calendar years.
    This is the quantity plotted in Fig. 8f.
    """
    sub = d[d.regime == regime] if regime else d
    tightm = (sub.rel.abs() < tight).astype(float)
    return tightm.rolling(win, center=True, min_periods=max(2, win // 2)).mean()

# ---------------- provenance / precision signature ---------------------------
def sigdigits(x):
    if not np.isfinite(x) or x==0: return np.nan
    s=('%.10g'%x).replace('-','').replace('.','').lstrip('0').rstrip('0')
    return len(s) if s else 1

def run(d, ref_mask):
    V=T1(d)
    par=fit_tolerance(d[ref_mask])
    Z,sigma=closure_Z(d,par)
    A=anomalies(d,ref_mask)
    T=temporal_score(d,ref_mask); S=spatial_score(d,A,ref_mask)
    src,SC=attribute(d,T,S,V)
    flag=pd.Series('N0-pass',index=d.index)
    flag[Z.abs()>2]='N1-suspect'
    flag[Z.abs()>3]='N2-fail'
    flag[V['any']]='N3-violation'
    return dict(V=V,par=par,Z=Z,sigma=sigma,A=A,T=T,S=S,src=src,SC=SC,flag=flag)

# ---------------- episode-level attribution ----------------------------------
def episodes(flag, min_level=('N1-suspect','N2-fail','N3-violation'), gap=1):
    bad = flag.isin(min_level).values
    idx = np.where(bad)[0]
    eps=[]
    if len(idx)==0: return eps
    start=idx[0]; prev=idx[0]
    for i in idx[1:]:
        if i-prev<=gap+1: prev=i
        else: eps.append((start,prev)); start=i; prev=i
    eps.append((start,prev))
    return eps

def attribute_episodes(d, SC, flag, Z=None):
    """Episode-level attribution.  `Z` (the standardised closure statistic) is
    optional; when given, maxZ is the largest |Z| inside each episode, which is
    the quantity reported in Table 4 of the manuscript."""
    if Z is None:
        Z = d['Z'] if 'Z' in d.columns else None
    out=[]
    for a,b in episodes(flag):
        sl=slice(a,b+1)
        tot=SC.iloc[sl].sum()
        share=(tot/tot.sum()).round(3)
        out.append(dict(start=d.index[a],end=d.index[b],n=b-a+1,
                        source=tot.idxmax(),
                        confidence=float(share.max()),
                        shares=share.to_dict(),
                        maxZ=(None if Z is None
                              else round(float(Z.iloc[sl].abs().max()),1))))
    return pd.DataFrame(out)
