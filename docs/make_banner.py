"""Draws docs/banner.png for the README from the processed record (run after step 01)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import ncqc
from paths import PROCESSED, ROOT

df = pd.read_csv(PROCESSED / 'merged.csv', parse_dates=['date']).set_index('date')
d = ncqc.build(df)
CLEAN = (d.yr <= 2024) & (~d.yr.isin([2016, 2018, 2019]))
r = ncqc.run(d, CLEAN)
flag = r['flag'] != 'N0-pass'

W, H = 16, 4.2
fig = plt.figure(figsize=(W, H), dpi=150)
bg = fig.add_axes([0, 0, 1, 1]); bg.axis('off')
grad = np.linspace(0, 1, 512)[None, :]
cmap = LinearSegmentedColormap.from_list('b', ['#06213d', '#0b3d63', '#11607e'])
bg.imshow(grad, aspect='auto', cmap=cmap, extent=[0, 1, 0, 1])
bg.set_xlim(0, 1); bg.set_ylim(0, 1)

# hydrographs on the right
ax = fig.add_axes([0.40, 0.30, 0.58, 0.62]); ax.axis('off')
cols = {'FS': '#e8f1f8', 'SP': '#f2a65a', 'LR': '#5fd3b4'}
for s, lw, a in [('FS', 2.2, 0.95), ('SP', 1.5, 0.9), ('LR', 1.5, 0.9)]:
    ax.plot(d.index, d['Q_' + s], color=cols[s], lw=lw, alpha=a, solid_capstyle='round')
ax.set_yscale('log'); ax.set_xlim(d.index[0], d.index[-1])
# closure-error strip with flagged months
ax2 = fig.add_axes([0.40, 0.10, 0.58, 0.16]); ax2.axis('off')
ax2.axhspan(-2, 2, color='white', alpha=0.10, lw=0)
ax2.plot(d.index, d.rel.clip(-25, 25), color='#cfe3f1', lw=1.0, alpha=0.9)
ax2.scatter(d.index[flag], d.rel.clip(-25, 25)[flag], s=10, color='#ff7a59', zorder=3, lw=0)
ax2.set_xlim(d.index[0], d.index[-1]); ax2.set_ylim(-27, 27)

# text block
bg.text(0.035, 0.66, 'NCQC', color='white', fontsize=58, fontweight='bold', family='DejaVu Sans', va='center')
bg.text(0.037, 0.43, 'Network-Constrained Quality Control\nfor river discharge records', color='#d7e7f3',
        fontsize=17, va='center', linespacing=1.3)
bg.text(0.037, 0.17, 'Mackenzie–Liard confluence, NWT, Canada  ·  2000–2026', color='#9fc3dc', fontsize=11.5, va='center')
# legend chips
for i, (lab, c) in enumerate([('Fort Simpson', cols['FS']), ('Strong Point', cols['SP']),
                              ('Liard', cols['LR']), ('flagged month', '#ff7a59')]):
    x = 0.41 + i * 0.105
    bg.plot([x, x + 0.018], [0.955, 0.955], color=c, lw=3, solid_capstyle='round')
    bg.text(x + 0.024, 0.955, lab, color='#d7e7f3', fontsize=9.5, va='center')
fig.savefig(ROOT / 'docs' / 'banner.png', dpi=150, facecolor=fig.get_facecolor())
print('banner ok')
