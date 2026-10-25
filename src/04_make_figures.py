"""
04_make_figures.py
==================
Generates Figures 1-9 of the manuscript plus summary.json and episodes.csv.

    Fig1  locator map: Mackenzie-Liard region, Canada, three-gauge confluence
    Fig2  hydrographs, annual regime, partition climatology
    Fig3  mass-balance closure: time series, scatter, error distribution
    Fig4  NDVI regime stratification, tolerance envelope, epoch non-stationarity
    Fig5  information content of candidate constraints
    Fig6  NCQC algorithm decision chain
    Fig7  standardized closure statistic, flag matrix, CUSUM, attribution shares
    Fig8  case studies (2016; 2025-26 provisional) + T0 and T6 diagnostics
    Fig9  detection and attribution skill from the injection experiment

Usage:  python src/04_make_figures.py   (needs merged.csv and injection.csv)
Output: figures/Fig1.png ... figures/Fig9.png (Fig6 also as PDF),
        data/processed/summary.json, data/processed/episodes.csv
"""

import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Polygon
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
from math import comb
import ncqc, os, warnings
from paths import PROCESSED, FIGURES
warnings.filterwarnings('ignore', module='mpl_toolkits.basemap')
from mpl_toolkits.basemap import Basemap

plt.rcParams.update({'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,
 'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,
 'figure.dpi':160,'savefig.dpi':300,'legend.frameon':False,'axes.grid':True,
 'grid.alpha':0.25,'grid.linewidth':0.5})
OUT=str(FIGURES)
C={'FS':'#1b3b6f','SP':'#c1440e','LR':'#1f8a70'}
CR={'Ice-affected':'#4a7fb5','Transition':'#e0a33c','Open-water':'#3f9b5a'}

df=pd.read_csv(PROCESSED / 'merged.csv',parse_dates=['date']).set_index('date')
d=ncqc.build(df)
CLEAN=(d.yr<=2024)&(~d.yr.isin([2016,2018,2019]))
r=ncqc.run(d,CLEAN)
Z=r['Z']; d['Z']=Z; d['flag']=r['flag']; d['src']=r['src']
E=ncqc.attribute_episodes(d,r['SC'],r['flag'],Z=Z)
Sp,Sn=ncqc.cusum(Z)
SHORT={'FS':'Fort Simpson (10GC001)','SP':'Strong Point (10FB006)',
       'LR':'Liard near the mouth (10ED002)'}
NAMES={'FS':'Mackenzie R. at Fort Simpson (10GC001)','SP':'Mackenzie R. at Strong Point (10FB006)',
       'LR':'Liard R. near the mouth (10ED002)'}


# ---------------------------------------------------------------------------
# schematic helpers (Fig. 1a planform, Fig. 5 decision chain)
# ---------------------------------------------------------------------------
WATER='#c3d8ea'

# --- Figure 1a locator map ---------------------------------------------------
# Black-and-white locator drawn from the boundary, coastline and river data
# bundled with basemap (GSHHS / WDB-II), so it needs no network access.
RED='#d62728'
STATIONS={'SP':(-120.7917,61.8165,'10FB006','Strong Point'),
          'LR':(-121.2278,61.7426,'10ED002','Liard near the mouth'),
          'FS':(-121.3589,61.8685,'10GC001','Fort Simpson')}
SITE=(-121.62,-120.55,61.705,61.905)      # confluence detail extent
REG =(-126.5,-111.5,59.6,64.2)            # Mackenzie-Liard region panel extent
LAT0=0.5*(SITE[2]+SITE[3])
KM_LAT=110.57; KM_LON=111.32*np.cos(np.radians(LAT0))
SIDE={'SP':1,'LR':1,'FS':-1}              # +1 label above the marker, -1 below
# locator palette: water/land carry the only colour, boundaries stay near-black
WATER='#cfe2f1'; LAND='#f7f5f0'; COAST='#3f5160'
BORDER='#16202b'; PROV='#9aa7b4'; RIVER='#4e8ab8'
HALO=[pe.withStroke(linewidth=2.1,foreground='w',alpha=.92)]

def _band(ax,pts,hw,fc=WATER,z=3):
    """Channel ribbon along a Bezier centreline; half-width hw (data units)."""
    P=np.array(pts,float); t=np.linspace(0,1,260); n=len(P)-1
    cur=sum(comb(n,i)*((1-t)**(n-i))[:,None]*(t**i)[:,None]*P[i] for i in range(n+1))
    d=np.gradient(cur,axis=0); d/=np.linalg.norm(d,axis=1)[:,None]
    nrm=np.column_stack([-d[:,1],d[:,0]])
    ax.add_patch(Polygon(np.vstack([cur+nrm*hw,(cur-nrm*hw)[::-1]]),closed=True,
                         fc=fc,ec='none',zorder=z))
    return cur,nrm

def _gauge(ax,cur,nrm,frac,hw,code,name,col,side=1):
    """Gauging-section marker across the channel, with station labels."""
    i=int(frac*(len(cur)-1)); x,y=cur[i]; nx,ny=nrm[i]; L=hw*1.55
    ax.plot([x-nx*L,x+nx*L],[y-ny*L,y+ny*L],color=col,lw=1.6,
            solid_capstyle='round',zorder=6)
    ax.plot([x],[y],marker='o',ms=4.4,mfc='w',mec=col,mew=1.5,zorder=7)
    o1=side*(L+0.30); o2=side*(L+0.12)
    ax.text(x+nx*o1,y+ny*o1,code,ha='center',va='center',fontsize=8.1,
            weight='bold',color=col,zorder=7)
    ax.text(x+nx*o2,y+ny*o2,name,ha='center',va='center',fontsize=6.8,
            color='0.38',zorder=7)

def _paint_site(ax):
    """Fig. 1c: the three-gauge confluence, gauges in red."""
    m=Basemap(projection='cyl',llcrnrlon=SITE[0],urcrnrlon=SITE[1],
              llcrnrlat=SITE[2],urcrnrlat=SITE[3],resolution='f',ax=ax,fix_aspect=False)
    m.drawmapboundary(fill_color=WATER,linewidth=1.0)
    m.fillcontinents(color=LAND,lake_color=WATER)
    m.drawcoastlines(linewidth=0.6,color=COAST)
    m.drawrivers(linewidth=0.9,color=RIVER)
    ax.set_aspect(KM_LAT/KM_LON); ax.grid(False)
    m.drawparallels([61.75,61.80,61.85,61.90],labels=[1,0,0,0],fontsize=5.8,
                    linewidth=0.25,color='0.85',dashes=[3,3])
    m.drawmeridians([-121.4,-121.2,-121.0,-120.8,-120.6],labels=[0,0,0,1],fontsize=5.8,
                    linewidth=0.25,color='0.85',dashes=[3,3])
    HALO=[pe.withStroke(linewidth=2.1,foreground='w',alpha=.92)]
    for k,(lon,lat,code,name) in STATIONS.items():
        ax.plot([lon],[lat],marker='o',ms=6.2,mfc=RED,mec='k',mew=0.8,zorder=10)
        dy=SIDE[k]*0.0105; va='bottom' if dy>0 else 'top'
        ax.text(lon,lat+dy,code,ha='center',va=va,fontsize=7.4,weight='bold',color='k',
                zorder=11,path_effects=HALO)
        ax.text(lon,lat+dy+SIDE[k]*0.0165,name,ha='center',va=va,fontsize=6.2,color='0.3',
                zorder=11,path_effects=HALO)
    dl=10.0/KM_LON; bx,by=SITE[0]+0.030,SITE[2]+0.017
    ax.add_patch(Rectangle((bx,by),dl,0.0050,fc='k',ec='k',lw=0.5,zorder=10))
    ax.add_patch(Rectangle((bx,by),dl/2,0.0050,fc='w',ec='k',lw=0.5,zorder=10))
    ax.text(bx+dl/2,by+0.0095,'10 km',ha='center',va='bottom',fontsize=5.9,zorder=10,
            path_effects=HALO)
    nx=SITE[1]-0.035
    ax.annotate('',xy=(nx,SITE[3]-0.024),xytext=(nx,SITE[3]-0.060),
                arrowprops=dict(arrowstyle='-|>',lw=0.9,color='k'),zorder=10)
    ax.text(nx,SITE[3]-0.020,'N',ha='center',va='bottom',fontsize=6.8,
            weight='bold',zorder=10,path_effects=HALO)


def _paint_region(ax):
    """Fig. 1a: the Mackenzie-Liard region, with the confluence boxed."""
    mr=Basemap(projection='cyl',llcrnrlon=REG[0],urcrnrlon=REG[1],
               llcrnrlat=REG[2],urcrnrlat=REG[3],resolution='h',ax=ax,fix_aspect=False)
    mr.drawmapboundary(fill_color=WATER,linewidth=0.9)
    mr.fillcontinents(color=LAND,lake_color=WATER)
    mr.drawcoastlines(linewidth=0.45,color=COAST)
    mr.drawrivers(linewidth=0.60,color=RIVER)
    mr.drawstates(linewidth=0.65,color=BORDER)
    ax.set_aspect(KM_LAT/(111.32*np.cos(np.radians(0.5*(REG[2]+REG[3]))))); ax.grid(False)
    ax.add_patch(Rectangle((SITE[0]-0.33,SITE[2]-0.24),(SITE[1]-SITE[0])+0.66,
                 (SITE[3]-SITE[2])+0.48,fc='none',ec=RED,lw=1.4,zorder=9))
    ax.text(0.035,0.035,'Mackenzie\u2013Liard region',transform=ax.transAxes,fontsize=6.6,
            color=BORDER,va='bottom',ha='left',zorder=11,
            bbox=dict(boxstyle='square,pad=0.2',fc='w',ec='none',alpha=.92))


def _paint_canada(ax):
    """Fig. 1b: Canada, with the study area boxed."""
    mc=Basemap(projection='lcc',lon_0=-96,lat_0=59,lat_1=49,lat_2=77,
               width=7.20e6,height=4.66e6,resolution='l',ax=ax)
    mc.drawmapboundary(fill_color=WATER,linewidth=0.9)
    mc.fillcontinents(color=LAND,lake_color=WATER)
    mc.drawcoastlines(linewidth=0.40,color=COAST)
    mc.drawstates(linewidth=0.40,color=PROV)
    mc.drawcountries(linewidth=0.85,color=BORDER)
    ax.grid(False)
    xr,yr=mc(-121.2,61.8)
    ax.plot(xr,yr,marker='s',ms=6.5,mfc='none',mec=RED,mew=1.6,zorder=9)
    ax.plot(xr,yr,marker='.',ms=2.4,color=RED,zorder=10)
    ax.text(0.035,0.035,'Canada',transform=ax.transAxes,fontsize=6.6,color=BORDER,
            va='bottom',ha='left',zorder=11,
            bbox=dict(boxstyle='square,pad=0.2',fc='w',ec='none',alpha=.92))


# ============================ FIGURE 1 : locator ============================
# The three panels are tiled flush: axes positions are computed so that each box
# aspect equals its map's own aspect, which is what lets them share edges without
# matplotlib shrinking an axes to satisfy set_aspect.
_FW,_LM,_RM,_BM,_TM=7.4,0.50,0.06,0.46,0.06      # inches: width and margins
_BW=_FW-_LM-_RM                                   # block width
_ASP_TOP=7.20/4.66                                # region and Canada panels
_ASP_BOT=(SITE[1]-SITE[0])*KM_LON/((SITE[3]-SITE[2])*KM_LAT)
_HT=(_BW/2)/_ASP_TOP                              # top-row height
_HB=_BW/_ASP_BOT                                  # bottom-row height
_FH=_BM+_HT+_HB+_TM
fig=plt.figure(figsize=(_FW,_FH))
_x0,_y0=_LM/_FW,_BM/_FH
_w2,_hT,_hB=(_BW/2)/_FW,_HT/_FH,_HB/_FH
_axC=fig.add_axes([_x0,_y0,_BW/_FW,_hB])          # (c) confluence, bottom
_axA=fig.add_axes([_x0,_y0+_hB,_w2,_hT])          # (a) region, top left
_axB=fig.add_axes([_x0+_w2,_y0+_hB,_w2,_hT])      # (b) Canada, top right
_paint_region(_axA); _paint_canada(_axB); _paint_site(_axC)
for _ax,_lab in [(_axA,'(a)'),(_axB,'(b)'),(_axC,'(c)')]:
    for _sp in _ax.spines.values():
        _sp.set_visible(True); _sp.set_linewidth(1.0); _sp.set_color('#16202b')
    _ax.text(0.012,0.985,_lab,transform=_ax.transAxes,ha='left',va='top',
             fontsize=9,weight='bold',zorder=20,
             bbox=dict(boxstyle='square,pad=0.22',fc='w',ec='#16202b',lw=0.7))
for _ax in (_axA,_axB):
    _ax.set_xticks([]); _ax.set_yticks([])
fig.savefig(f'{OUT}/Fig1.png',dpi=300,facecolor='w'); plt.close(fig)
print('fig1 ok')

# ============================ FIGURE 2 : hydrology ==========================
fig=plt.figure(figsize=(7.4,5.1))
gs=fig.add_gridspec(2,2,height_ratios=[1.15,1],hspace=0.42,wspace=0.28)
ax=fig.add_subplot(gs[0,:])
for s in ['FS','SP','LR']:
    ax.plot(d.index,d['Q_'+s],lw=.85,color=C[s],label=SHORT[s])
ax.set_ylabel('Monthly $Q$  (m$^3$ s$^{-1}$)'); ax.set_yscale('log')
ax.legend(fontsize=6.5,ncol=3,loc='lower left',bbox_to_anchor=(0.075,1.005),
          frameon=False,handlelength=1.6,columnspacing=1.6)   # above the axes: it
                                                              # overlapped the peaks
ax.set_title('(a)',loc='left')
ax.axvspan(pd.Timestamp('2025-01-01'),d.index[-1],color='0.85',alpha=.5,zorder=0)
ax.text(pd.Timestamp('2025-03-01'),1.6e4,'provisional',fontsize=6.5,color='0.35')

ax=fig.add_subplot(gs[1,0])
cl=d.groupby('m')[['Q_FS','Q_SP','Q_LR']].median()
ax.fill_between(cl.index,0,cl.Q_SP,color=C['SP'],alpha=.55,label='Strong Point')
ax.fill_between(cl.index,cl.Q_SP,cl.Q_SP+cl.Q_LR,color=C['LR'],alpha=.55,label='Liard')
ax.plot(cl.index,cl.Q_FS,'k-',lw=1.5,label='Fort Simpson (obs.)')
ax.set_xticks(range(1,13)); ax.set_xticklabels(list('JFMAMJJASOND'))
ax.set_ylabel('m$^3$ s$^{-1}$')
ax.set_ylim(0,cl.Q_FS.max()*1.42)   # headroom for the legend
ax.legend(fontsize=6.4,loc='upper left')
ax.set_title('(b)',loc='left')

ax=fig.add_subplot(gs[1,1])
f=(d.Q_LR/(d.Q_SP+d.Q_LR))
bp=ax.boxplot([f[d.m==m].values for m in range(1,13)],widths=.6,showfliers=False,
              patch_artist=True,medianprops=dict(color='k',lw=1))
for p in bp['boxes']: p.set_facecolor(C['LR']); p.set_alpha(.4); p.set_edgecolor('0.4')
ax.set_xticks(range(1,13)); ax.set_xticklabels(list('JFMAMJJASOND'))
ax.set_ylabel('Liard share $f_t$  (–)')
ax.set_title('(c)',loc='left')
fig.savefig(f'{OUT}/Fig2.png',bbox_inches='tight'); plt.close(fig)
print('fig2 ok')

# ============================ FIGURE 3 =======================================
fig=plt.figure(figsize=(7.4,4.6))
gs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=0.55,wspace=0.28)
axs=[fig.add_subplot(gs[1,0]),fig.add_subplot(gs[0,:]),fig.add_subplot(gs[1,1])]
ax=axs[0]
for rg in CR:
    g=d[d.regime==rg]
    ax.scatter(g.Qsum,g.Q_FS,s=8,color=CR[rg],alpha=.7,lw=0,label=rg)
lim=[1e3,2.6e4]; ax.plot(lim,lim,'k--',lw=.8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('$Q_{SP}+Q_{LR}$  (m$^3$ s$^{-1}$)'); ax.set_ylabel('$Q_{FS}$  (m$^3$ s$^{-1}$)')
ax.legend(fontsize=6.2,loc='upper left'); ax.set_title('(b)',loc='left')
ax.text(.97,.06,'$r^2$ = %.4f\nNSE = %.4f'%(np.corrcoef(d.Qsum,d.Q_FS)[0,1]**2,
        1-((d.Q_FS-d.Qsum)**2).sum()/((d.Q_FS-d.Q_FS.mean())**2).sum()),
        transform=ax.transAxes,ha='right',fontsize=6.8)

ax=axs[1]
# shade NDVI-classified ice-affected months (T2) so that winter closure can be read directly
_ice=(d.regime=='Ice-affected').values
for _t,_f in zip(d.index,_ice):
    if _f: ax.axvspan(_t,_t+pd.offsets.MonthBegin(1),color=CR['Ice-affected'],alpha=.16,lw=0,zorder=0)
ax.axhspan(-2,2,color='0.85',alpha=.6,lw=0)
ax.plot(d.index,d.rel,lw=.8,color='0.35')
sel=d.flag!='N0-pass'
ax.scatter(d.index[sel],d.rel[sel],s=12,color='#c1440e',zorder=3,lw=0)
ax.axhline(0,color='k',lw=.6)
ax.set_ylabel('closure error $r_t/Q_{FS}$  (%)')
ax.set_title('(a)',loc='left')
ax.set_ylim(-45,26)
ax.legend(handles=[Rectangle((0,0),1,1,color=CR['Ice-affected'],alpha=.3,lw=0),
                   Rectangle((0,0),1,1,color='0.85',lw=0),
                   Line2D([],[],marker='o',ls='',ms=3.5,color='#c1440e')],
          labels=['ice-affected month (T2)','$\\pm$2 %','NCQC-flagged month'],
          fontsize=6.2,loc='lower right',bbox_to_anchor=(1.0,0.99),ncol=3,borderaxespad=0.)
ax.annotate('2016 episode',xy=(pd.Timestamp('2016-08-01'),17.4),xytext=(pd.Timestamp('2010-06-01'),19),
            fontsize=6.8,arrowprops=dict(arrowstyle='->',lw=.7))
ax.annotate('provisional ice-season bias',xy=(pd.Timestamp('2025-01-01'),-38.9),
            xytext=(pd.Timestamp('2014-01-01'),-40),fontsize=6.8,
            arrowprops=dict(arrowstyle='->',lw=.7))

ax=axs[2]
x=d.rel[CLEAN]
ax.hist(x,bins=np.arange(-8,8.5,0.5),color='0.6',edgecolor='w',lw=.4,density=True)
xs=np.linspace(-8,8,300)
sd=ncqc.mad(x)
ax.plot(xs,1/(sd*np.sqrt(2*np.pi))*np.exp(-xs**2/(2*sd**2)),'k-',lw=1.2,
        label='Gaussian, $\\sigma_{MAD}$=%.2f %%'%sd)
ax.set_xlabel('closure error (%)'); ax.set_ylabel('density')
ax.legend(fontsize=6.2); ax.set_title('(c)',loc='left')
fig.savefig(f'{OUT}/Fig3.png',bbox_inches='tight'); plt.close(fig)
print('fig2-3 ok')

# ============================ FIGURE 4 =======================================
fig=plt.figure(figsize=(7.4,4.8))
gs=fig.add_gridspec(2,3,hspace=0.55,wspace=0.35)
ax=fig.add_subplot(gs[0,:2])
mo=d.groupby('m')[['ndvi_FS','ndvi_SP','ndvi_LR']].median()
for s in ['FS','SP','LR']:
    ax.plot(mo.index,mo['ndvi_'+s],'o-',ms=3,lw=1,color=C[s],label=s)
ax.axhline(0.15,color='0.4',ls=':',lw=.9)
ax.text(1.1,0.17,'NDVI = 0.15',fontsize=6.4,color='0.35')
rgmo=d.groupby('m')['regime'].agg(lambda x:x.mode()[0])
for m,rg in rgmo.items():
    ax.axvspan(m-.5,m+.5,color=CR[rg],alpha=.12,lw=0,zorder=0)
ax.set_xticks(range(1,13)); ax.set_xticklabels(list('JFMAMJJASOND')); ax.set_xlim(.5,12.5)
ax.set_ylabel('MODIS NDVI at gauge pixel'); ax.legend(fontsize=6.4,ncol=3,loc='upper left')
ax.set_title('(a)',loc='left')

ax=fig.add_subplot(gs[0,2])
for rg in CR:
    g=d[(d.regime==rg)&CLEAN]
    ax.scatter(g.Q_FS,g.res.abs(),s=7,color=CR[rg],alpha=.7,lw=0)
    q=np.array([g.Q_FS.min(),g.Q_FS.max()])
    k=r['par'][rg]['k']; s0=r['par'][rg]['s0']
    ax.plot(q,2*(k*q+s0),color=CR[rg],lw=1.2)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('$Q_{FS}$  (m$^3$ s$^{-1}$)'); ax.set_ylabel('$|r_t|$  (m$^3$ s$^{-1}$)')
ax.legend(handles=[Line2D([],[],color=CR[k],lw=1.4,marker='o',ms=3,ls='-',label=k)
                   for k in ['Ice-affected','Transition','Open-water']],
          fontsize=5.6,loc='upper left',handlelength=1.4)
ax.set_title('(b)',loc='left')

ax=fig.add_subplot(gs[1,:])
ep=pd.cut(d.yr,[1999,2009,2017,2019,2024,2026],labels=['2000–09','2010–17','2018–19','2020–24','2025–26'])
d['epoch']=ep
pos=0; ticks=[]; lab=[]
for e in ['2000–09','2010–17','2018–19','2020–24','2025–26']:
    for rg in ['Ice-affected','Transition','Open-water']:
        v=d.loc[(d.epoch==e)&(d.regime==rg),'rel'].dropna().values
        if len(v)<3: pos+=1; continue
        bp=ax.boxplot([v],positions=[pos],widths=.72,patch_artist=True,showfliers=True,
                      flierprops=dict(ms=2,mfc='0.4',mec='none'),medianprops=dict(color='k',lw=1.1))
        bp['boxes'][0].set_facecolor(CR[rg]); bp['boxes'][0].set_alpha(.55); bp['boxes'][0].set_edgecolor('0.4')
        pos+=1
    ticks.append(pos-2); lab.append(e); pos+=1
ax.set_xticks(ticks); ax.set_xticklabels(lab); ax.axhline(0,color='k',lw=.6)
ax.set_ylabel('closure error (%)'); ax.set_ylim(-42,22)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=CR[k],alpha=.55,label=k) for k in ['Ice-affected','Transition','Open-water']],
          fontsize=6.4,ncol=3,loc='lower left')
ax.set_title('(c)',loc='left',fontsize=8.5)
fig.savefig(f'{OUT}/Fig4.png',bbox_inches='tight'); plt.close(fig)

# ============================ FIGURE 5 =======================================
A=r['A']
fig,axs=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw=dict(wspace=0.52))
ax=axs[0]
ax.scatter(A.a_SP,A.a_FS,s=8,color=C['SP'],alpha=.6,lw=0,label='$a_{SP}$ (Strong Point)')
ax.scatter(A.a_LR,A.a_FS,s=8,color=C['LR'],alpha=.6,lw=0,label='$a_{LR}$ (Liard)')
for s,c in [('SP',C['SP']),('LR',C['LR'])]:
    b=np.polyfit(A['a_'+s],A.a_FS,1); xs=np.linspace(A['a_'+s].min(),A['a_'+s].max(),20)
    ax.plot(xs,np.polyval(b,xs),color=c,lw=1.1)
ax.set_xlabel("tributary log-anomaly  $a_{SP}$,  $a_{LR}$")
ax.set_ylabel("$a_{FS}$")
ax.legend(fontsize=6.2,loc='lower right'); ax.set_title('(a)',loc='left')
ax.text(.03,.95,'$r$ = %.2f (SP), %.2f (LR)'%(A.a_FS.corr(A.a_SP),A.a_FS.corr(A.a_LR)),
        transform=ax.transAxes,fontsize=6.6)

ax=axs[1]
for s in ['FS','SP','LR']:
    ar=d['runoff_'+s]-d.groupby('m')['runoff_'+s].transform('median')
    ax.scatter(ar,A['a_'+s],s=7,color=C[s],alpha=.55,lw=0,label='%s ($r$=%.2f)'%(s,ar.corr(A['a_'+s])))
_arr=np.concatenate([(d['runoff_'+s]-d.groupby('m')['runoff_'+s].transform('median')).dropna().values
                     for s in ['FS','SP','LR']])
_lo,_hi=np.percentile(_arr,[0.5,99.0]); _pad=0.08*(_hi-_lo)
ax.set_xlim(_lo-_pad,_hi+_pad)      # a few extreme months would otherwise flatten the cloud
ax.set_xlabel('ERA5-Land runoff anomaly (m)'); ax.set_ylabel('$Q$ log-anomaly')
ax.legend(fontsize=6.2,loc='lower right'); ax.set_title('(b)',loc='left')

ax=axs[2]
pred=['Mass balance\n($Q_{SP}+Q_{LR}$)','Neighbour\nanomaly (SP)','Neighbour\nanomaly (LR)',
      'Temporal\npersistence','ERA5-Land\nrunoff','ERA5-Land\nprecip.','MODIS\nNDVI']
lq=np.log(d.Q_FS)
def r2(x,y):
    m=np.isfinite(x)&np.isfinite(y)
    return np.corrcoef(x[m],y[m])[0,1]**2
vals=[r2(np.log(d.Qsum),lq),
      r2(A.a_SP,A.a_FS), r2(A.a_LR,A.a_FS),
      r2(np.log(d.Q_FS.shift(1)),lq),
      r2((d.runoff_FS-d.groupby('m')['runoff_FS'].transform('median')).values,A.a_FS.values),
      r2((d.precip_FS-d.groupby('m')['precip_FS'].transform('median')).values,A.a_FS.values),
      r2(d.ndvi_FS.values,A.a_FS.values)]
cols=['#1b3b6f','#3f6fa8','#3f6fa8','#7d8ca3','#b0b7c2','#b0b7c2','#b0b7c2']
_o=np.argsort(vals)[::-1]                      # rank by value so the bars read as a ranking
vals=[vals[i] for i in _o]; pred=[pred[i] for i in _o]; cols=[cols[i] for i in _o]
ax.barh(range(7),vals,color=cols); ax.set_yticks(range(7)); ax.set_yticklabels(pred,fontsize=6.2)
ax.invert_yaxis(); ax.set_xlabel('variance of $Q_{FS}$ (anomaly) explained, $R^2$')
for i,v in enumerate(vals):                    # never round 0.995 up to a misleading 1.00
    ax.text(v+.01,i,('%.2f'%v) if v<0.995 else ('%.3f'%v),va='center',fontsize=6.2)
ax.set_xlim(0,1.18); ax.set_title('(c)',loc='left')
fig.savefig(f'{OUT}/Fig5.png',bbox_inches='tight'); plt.close(fig)
print('fig4-5 ok')

# ============================ FIGURE 6 : decision chain ======================
# Drawn in a TikZ-like idiom: rectangles for processes, diamonds for decisions,
# stadia for terminals, orthogonal arrows, serif type.  Coordinates are inches.
_prev_family=plt.rcParams['font.family']
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral','DejaVu Serif'],
                     'mathtext.fontset':'stix'})
INK='#000000'; LW=1.0
FILL_P='#eceff4'; FILL_D='#f5eee1'; FILL_T='#ffffff'
FW,FH=6.4,8.6
fig=plt.figure(figsize=(FW,FH))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,FW); ax.set_ylim(0,FH)
ax.axis('off'); ax.set_aspect('equal')

def process(x,y,w,h,head,sub):
    ax.add_patch(Rectangle((x-w/2,y-h/2),w,h,fc=FILL_P,ec=INK,lw=LW,zorder=3,joinstyle='miter'))
    ax.text(x,y+0.10,head,ha='center',va='center',fontsize=11.0,color=INK,zorder=4)
    ax.text(x,y-0.12,sub,ha='center',va='center',fontsize=9.6,color=INK,zorder=4)

def decision(x,y,w,h,txt):
    ax.add_patch(Polygon([[x,y+h/2],[x+w/2,y],[x,y-h/2],[x-w/2,y]],closed=True,
        fc=FILL_D,ec=INK,lw=LW,zorder=3,joinstyle='miter'))
    ax.text(x,y,txt,ha='center',va='center',fontsize=11.0,color=INK,zorder=4)

def terminal(x,y,w,h,txt,fs=11.0):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,
        boxstyle='round,pad=0,rounding_size=%.3f'%(h/2),fc=FILL_T,ec=INK,lw=LW,zorder=3))
    ax.text(x,y,txt,ha='center',va='center',fontsize=fs,color=INK,zorder=4)

def arrow(p0,p1):
    ax.add_patch(FancyArrowPatch(p0,p1,arrowstyle='-|>,head_length=4.2,head_width=2.0',lw=LW,
                 color=INK,shrinkA=0,shrinkB=0,zorder=2,joinstyle='miter',capstyle='butt'))

CX,W=2.62,3.62; XF,WF,HF=5.62,0.74,0.34
Y={'in':8.22,'t0':7.52,'t1':6.60,'t2':5.70,'t3':4.78,'t4':3.74,'t5':2.98,'t6':2.00,'out':1.00}
HP,HD,HT=0.46,0.72,0.34

terminal(CX,Y['in'],W,0.40,'inputs:  $Q$ at the gauges,  NDVI,  ERA5-Land',fs=10.0)
process (CX,Y['t0'],W,HP,'T0   provenance signature','significant-digit count of reported values')
decision(CX,Y['t1'],W,HD,'T1   admissible?')
process (CX,Y['t2'],W,HP,'T2   regime classification','ice-affected / transition / open-water')
decision(CX,Y['t3'],W,HD,r'T3   $|Z_t| \leq 2$?')
process (CX,Y['t4'],W,HP,'T4   attribution','temporal and spatial scores, by station')
process (CX,Y['t5'],W,HP,'T5   episode and change point',r'contiguous flags;  CUSUM on $Z_t$')
decision(CX,Y['t6'],W,HD,'T6   closure independent?')
terminal(CX,Y['out'],W,0.40,'output:  flag,  attributed source,  confidence',fs=10.0)

for y,lab in [(Y['t1'],'N3'),(Y['t3'],'N0'),(Y['t5'],'N4'),(Y['t6'],'N5')]:
    terminal(XF,y,WF,HF,lab)
terminal(0.62,4.26,0.98,HF,'N1 / N2',fs=10.0)

# ---- spine -----------------------------------------------------------------
spine=[(Y['in']-0.20,Y['t0']+HP/2),(Y['t0']-HP/2,Y['t1']+HD/2),
       (Y['t1']-HD/2,Y['t2']+HP/2),(Y['t2']-HP/2,Y['t3']+HD/2),
       (Y['t3']-HD/2,Y['t4']+HP/2),(Y['t4']-HP/2,Y['t5']+HP/2),
       (Y['t5']-HP/2,Y['t6']+HD/2),(Y['t6']-HD/2,Y['out']+0.20)]
for a,b in spine: arrow((CX,a),(CX,b))
for y,lab in [(Y['t1'],'yes'),(Y['t3'],'no'),(Y['t6'],'yes')]:
    ax.text(CX+0.09,y-HD/2-0.17,lab,fontsize=9.0,color=INK,ha='left',va='center',zorder=5)

# ---- branches to the flag terminals ----------------------------------------
for y,lab in [(Y['t1'],'no'),(Y['t3'],'yes'),(Y['t6'],'no')]:
    arrow((CX+W/2,y),(XF-WF/2,y))
    ax.text((CX+W/2+XF-WF/2)/2,y+0.07,lab,fontsize=9.0,color=INK,ha='center',
            va='bottom',zorder=5)
arrow((CX+W/2,Y['t5']),(XF-WF/2,Y['t5']))
ax.text((CX+W/2+XF-WF/2)/2,Y['t5']+0.07,'step / drift',fontsize=9.0,color=INK,
        ha='center',va='bottom',zorder=5)
# N1 / N2 taken off the |Z| > 2 path
ax.plot([CX,CX],[Y['t3']-HD/2,4.26],color=INK,lw=LW,zorder=2,solid_capstyle='butt')
arrow((CX,4.26),(0.62+0.49,4.26))


fig.savefig(f'{OUT}/Fig6.png',facecolor='w',bbox_inches='tight',pad_inches=0.05,dpi=600)
fig.savefig(f'{OUT}/Fig6.pdf',facecolor='w',bbox_inches='tight',pad_inches=0.05,
            metadata={'CreationDate': None})
# trim the empty margin left below the output terminal
from PIL import Image, ImageChops
_im = Image.open(f'{OUT}/Fig6.png').convert('RGB')
_bb = ImageChops.difference(_im, Image.new('RGB', _im.size, (255, 255, 255))).getbbox()
_im.crop((max(0, _bb[0]-60), max(0, _bb[1]-60), min(_im.size[0], _bb[2]+60),
          min(_im.size[1], _bb[3]+60))).save(f'{OUT}/Fig6.png', dpi=(600, 600))
plt.close(fig)
plt.rcParams['font.family']=_prev_family
plt.rcParams['mathtext.fontset']='dejavusans'
print('fig6 ok')

# ============================ FIGURE 7 =======================================
FCOL={'N0-pass':'#e9eef2','N1-suspect':'#f0c674','N2-fail':'#d9743c','N3-violation':'#a32b12'}
fig=plt.figure(figsize=(7.4,5.6))
gs=fig.add_gridspec(3,2,height_ratios=[1,1.15,1],hspace=0.62,wspace=0.26)
ax=fig.add_subplot(gs[0,:])
ax.axhspan(-2,2,color='#d9e4d2',alpha=.5,lw=0); ax.axhspan(-3,-2,color='#f7ecd2',alpha=.7,lw=0)
ax.axhspan(2,3,color='#f7ecd2',alpha=.7,lw=0)
ax.plot(d.index,Z,lw=.8,color='0.3')
for f in ['N1-suspect','N2-fail','N3-violation']:
    m=d.flag==f
    ax.scatter(d.index[m],Z[m],s=14,color=FCOL[f],ec='0.3',lw=.3,zorder=3,label=f)
ax.set_ylabel('$Z_t$  (–)'); ax.set_ylim(-30,14); ax.legend(fontsize=6.3,ncol=3,loc='lower left')
ax.set_title('(a)',loc='left')

ax=fig.add_subplot(gs[1,0])
piv=d.pivot_table(index='m',columns='yr',values='Z',aggfunc='first')
lev=pd.DataFrame(np.select([piv.abs()<=2,piv.abs()<=3,piv.abs()<=99],[0,1,2]),index=piv.index,columns=piv.columns)
lev=lev.where(piv.notna())
for t in d.index[r['V']['any'].values]: lev.loc[t.month,t.year]=3
from matplotlib.colors import ListedColormap, BoundaryNorm
cm=ListedColormap([FCOL['N0-pass'],FCOL['N1-suspect'],FCOL['N2-fail'],FCOL['N3-violation']])
ax.pcolormesh(piv.columns,piv.index,lev.values,cmap=cm,norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4),
              edgecolors='w',lw=.3)
ax.set_yticks(range(1,13)); ax.set_yticklabels(list('JFMAMJJASOND'),fontsize=6.5)
ax.invert_yaxis(); ax.set_title('(b)',loc='left'); ax.grid(False)

ax=fig.add_subplot(gs[1,1])
ax.plot(d.index,Sp,color='#c1440e',lw=1,label='CUSUM$^+$')
ax.plot(d.index,Sn,color='#1b3b6f',lw=1,label='CUSUM$^-$')
ax.axhline(10,ls=':',color='0.4',lw=.8); ax.axhline(-10,ls=':',color='0.4',lw=.8)
ax.text(d.index[4],-24,'decision interval  $h=\\pm10$',fontsize=5.8,color='0.4')
ax.set_ylabel('CUSUM of $Z_t$'); ax.legend(fontsize=6.3,loc='upper left')
ax.set_yticks([-75,-50,-25,0,25])   # CUSUM- bottoms out at -80.7, so label -25 and -75
ax.set_title('(c)',loc='left')
ax.annotate('2018–19 step',xy=(pd.Timestamp('2019-10-01'),41),xytext=(pd.Timestamp('2008-01-01'),34),
            fontsize=6.4,arrowprops=dict(arrowstyle='->',lw=.7))

ax=fig.add_subplot(gs[2,:])
Eb=E[E.n>=1].copy()
Eb['lab']=[f"{a:%Y-%m}" for a in Eb.start]
x=np.arange(len(Eb)); bot=np.zeros(len(Eb))
for s in ['FS','SP','LR']:
    v=np.array([sh[s] for sh in Eb.shares])
    ax.bar(x,v,bottom=bot,color=C[s],alpha=.8,label=NAMES[s].split('(')[0].strip(),width=.72)
    bot+=v
ax.set_xticks(x); ax.set_xticklabels(Eb.lab,rotation=90,fontsize=5.8)
ax.set_ylabel('attribution share'); ax.set_ylim(0,1)
for xi,c in zip(x,Eb.confidence):
    ax.text(xi,1.02,'%.2f'%c,ha='center',fontsize=5.2,color='0.35')
ax.legend(fontsize=6.2,ncol=3,loc='upper center',bbox_to_anchor=(.5,-0.62))
ax.set_title('(d)',loc='left')
fig.savefig(f'{OUT}/Fig7.png',bbox_inches='tight'); plt.close(fig)
print('fig7 ok')

# ============================ FIGURE 8 =======================================
fig=plt.figure(figsize=(7.4,6.2))
gs=fig.add_gridspec(3,2,hspace=0.62,wspace=0.28,height_ratios=[1,1,1])

def case(axQ,axR,t0,t1,title):
    sub=d.loc[t0:t1]
    axQ.plot(sub.index,sub.Q_FS,'-o',ms=2.5,lw=1.1,color=C['FS'],label='$Q_{FS}$ reported')
    axQ.plot(sub.index,sub.Qsum,'--s',ms=2.5,lw=1.1,color='0.35',label='$Q_{SP}+Q_{LR}$')
    axQ.plot(sub.index,sub.Q_SP,'-',lw=.9,color=C['SP'],alpha=.8,label='$Q_{SP}$')
    axQ.plot(sub.index,sub.Q_LR,'-',lw=.9,color=C['LR'],alpha=.8,label='$Q_{LR}$')
    axQ.set_yscale('log'); axQ.set_ylabel('m$^3$ s$^{-1}$')
    _lo=float(sub[['Q_FS','Q_SP','Q_LR']].min().min())
    axQ.set_ylim(bottom=_lo/3.2)          # clear band under the data for the legend
    axQ.legend(fontsize=5.9,ncol=2,loc='lower left')
    axQ.set_title(title,loc='left')
    cols=[CR[x] for x in sub.regime]
    axR.bar(sub.index,sub.rel,width=22,color=cols,edgecolor='0.35',lw=.4)
    axR.axhline(0,color='k',lw=.6)
    for t in sub.index[ncqc.T1(sub)['any'].values]:
        axR.plot(t,sub.rel.loc[t]*1.06,marker='v',ms=4,color='#a32b12')
    axR.set_ylabel('closure error (%)')
    import matplotlib.dates as mdates
    for a in (axQ,axR):
        a.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
        a.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        for lb in a.get_xticklabels(): lb.set_fontsize(6.0); lb.set_rotation(45); lb.set_ha('right')

a1=fig.add_subplot(gs[0,0]); a2=fig.add_subplot(gs[1,0])
case(a1,a2,'2015-06','2017-06','(a)')
a2.set_title('(c)',loc='left',fontsize=8)
a3=fig.add_subplot(gs[0,1]); a4=fig.add_subplot(gs[1,1])
case(a3,a4,'2023-07','2026-07','(b)')
a4.set_title('(d)',loc='left',fontsize=8)

ax=fig.add_subplot(gs[2,0])
sd=pd.DataFrame({s:d['Q_'+s].apply(ncqc.sigdigits) for s in ['FS','SP','LR']},index=d.index)
for s in ['FS','SP','LR']:
    ax.plot(d.index,sd[s].rolling(6,min_periods=3).mean(),lw=1.1,color=C[s],label=s)
ax.axhline(3,ls=':',color='0.4',lw=.8)
ax.text(d.index[10],3.25,'3 significant digits (approved / rounded)',fontsize=6.2,color='0.35')
ax.set_ylabel('significant digits\n(6-month mean)'); ax.legend(fontsize=6.3,ncol=3)
ax.set_title('(e)',loc='left')

ax=fig.add_subplot(gs[2,1])
ice=d[d.regime=='Ice-affected']
tight=(ice.rel.abs()<0.5).astype(float).rolling(12,min_periods=6,center=True).mean()
ax.fill_between(tight.index,0,tight.values,color='#7a5fa6',alpha=.35,lw=0)
ax.plot(tight.index,tight.values,color='#7a5fa6',lw=1.1)
ax.axhline(0.7,ls='--',color='0.3',lw=.8)
ax.text(d.index[6],0.74,'N5 threshold',fontsize=6.2,color='0.3')
ax.set_ylabel('ice-season months with\n$|r_t/Q_{FS}|<0.5$ %  (fraction)'); ax.set_ylim(0,1.05)
ax.set_title('(f)',loc='left')
fig.savefig(f'{OUT}/Fig8.png',bbox_inches='tight'); plt.close(fig)
print('fig8 ok')

# ============================ FIGURE 9 =======================================
R=pd.read_csv(PROCESSED / 'injection.csv')
fig=plt.figure(figsize=(7.4,4.8))
gs8=fig.add_gridspec(2,2,hspace=0.55,wspace=0.32)
axs=[fig.add_subplot(gs8[0,0]),fig.add_subplot(gs8[0,1]),fig.add_subplot(gs8[1,0]),fig.add_subplot(gs8[1,1])]
TC={'bias':'#1b3b6f','spike':'#c1440e','step':'#1f8a70','frozen':'#7a5fa6'}
ax=axs[0]
for typ,g in R.groupby('typ'):
    p=g.groupby('delta')['hit'].mean()
    ax.plot(p.index*100,p.values,'-o',ms=3,lw=1.1,color=TC[typ],label=typ)
ax.axhline(0.5,ls=':',color='0.5',lw=.8)
ax.set_xlabel('injected error $\\delta$ (%)'); ax.set_ylabel('probability of detection')
ax.set_ylim(0,1.05); ax.legend(fontsize=6.2,loc='lower right')
ax.text(0.015,0.045,'$\\delta$ is nominal for frozen-value trials',
        transform=ax.transAxes,fontsize=5.7,color='0.4',style='italic')
ax.set_title('(a)',loc='left')

ax=axs[1]
for rg,g in R[R.typ!='frozen'].groupby('regime'):
    p=g.groupby('delta')['hit'].mean()
    ax.plot(p.index*100,p.values,'-o',ms=3,lw=1.1,color=CR[rg],label=rg)
ax.axhline(0.5,ls=':',color='0.5',lw=.8)
ax.set_xlabel('injected error $\\delta$ (%)'); ax.set_ylabel('probability of detection')
ax.set_ylim(0,1.05); ax.legend(fontsize=6.0,loc='lower right')
ax.set_title('(b)',loc='left')

ax=axs[2]
levs=np.arange(1.0,5.01,0.1)
clean_rate=[]
base=pd.read_csv(PROCESSED / 'merged.csv',parse_dates=['date']).set_index('date')
db=ncqc.build(base); Zb,_=ncqc.closure_Z(db,r['par']); Vb=ncqc.T1(db)
for lev in levs:
    clean_rate.append(float(((Zb.abs()>lev)|Vb['any'])[CLEAN].mean()))
for dl,col in [(0.05,'#9ab5d6'),(0.10,'#4f7fb5'),(0.20,'#1b3b6f')]:
    g=R[np.isclose(R.delta,dl)]
    pod=[float(((g.maxZ>lev)|g.viol).mean()) for lev in levs]
    ax.plot(clean_rate,pod,'-',lw=1.3,color=col,label='$\\delta$ = %d %%'%(dl*100))
for dl,col in [(0.05,'#9ab5d6'),(0.10,'#4f7fb5'),(0.20,'#1b3b6f')]:
    g=R[np.isclose(R.delta,dl)]; i=np.argmin(np.abs(levs-2.0))
    ax.plot(clean_rate[i],float(((g.maxZ>levs[i])|g.viol).mean()),'o',ms=5,mfc='w',mec=col,mew=1.3)
ax.plot([0,.3],[0,.3],'k--',lw=.7)
ax.text(0.045,0.055,'chance ($y=x$)',fontsize=5.8,color='0.35',ha='left',va='bottom',
        rotation=11,rotation_mode='anchor')
ax.set_xlabel('false-alarm rate on clean record'); ax.set_ylabel('probability of detection')
ax.set_xlim(0,0.30); ax.set_ylim(0,1.05)
ax.legend(fontsize=6.2,loc='lower right',bbox_to_anchor=(1.0,0.31))   # above the diagonal
ax.set_title('(c)',loc='left')

ax=axs[3]
acc=R[R.hit].groupby('station')['correct'].mean().reindex(['SP','LR','FS'])
n=R[R.hit].groupby('station').size().reindex(['SP','LR','FS'])
ax.barh(range(3),acc.values,color=[C[s] for s in acc.index],alpha=.85)
ax.axvline(1/3,ls='--',color='k',lw=.8)
ax.set_yticks(range(3))
ax.set_yticklabels(['SP\nStrong Point\n(main stem)','LR\nLiard\n(tributary)',
                    'FS\nFort Simpson\n(closure)'],fontsize=5.8,linespacing=1.25)
ax.invert_yaxis(); ax.set_xlim(0,1.0); ax.set_xlabel('correct attribution rate')
for i,(v,k) in enumerate(zip(acc.values,n.values)):
    ax.text(v+.02,i,'%.2f  (n = %d)'%(v,k),va='center',fontsize=6.4)
ax.text(1/3+.015,-0.44,'chance',fontsize=6.2,color='0.35')
ax.set_title('(d)',loc='left')
fig.savefig(f'{OUT}/Fig9.png',bbox_inches='tight'); plt.close(fig)
print('fig9 ok')

# ---- summary numbers for the manuscript ------------------------------------
summ={}
summ['n_months']=len(d)
summ['flag_counts']=d.flag.value_counts().to_dict()
summ['par']=r['par']
summ['pod']=R.groupby(['typ','delta'])['hit'].mean().round(3).to_dict()
summ['pod_regime']=R[R.typ!='frozen'].groupby(['regime','delta'])['hit'].mean().round(3).to_dict()
summ['att']=R[R.hit].groupby('station')['correct'].mean().round(3).to_dict()
summ['far']={float(l):float(((Zb.abs()>l)|Vb['any'])[CLEAN].mean()) for l in [1.5,2,2.5,3,4]}
import json; open(PROCESSED / 'summary.json','w').write(json.dumps({k:{str(kk):vv for kk,vv in v.items()} if isinstance(v,dict) else v for k,v in summ.items()},indent=1,default=str))
E.to_csv(PROCESSED / 'episodes.csv',index=False)
print(E.to_string())
