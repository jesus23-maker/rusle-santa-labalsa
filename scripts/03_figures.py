import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from PIL import Image
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8})
SRC='data/layers'; OUT='results'
m=json.load(open(f'{SRC}/meta.json')); W,S,E,N=m['bounds']; ext=[W,E,S,N]
basin=np.array(m['basin']['coordinates'][0])
lat0=(S+N)/2; aspect=1/np.cos(np.radians(lat0))
def img(n): return np.array(Image.open(f'{SRC}/{n}.png').convert('RGBA'))
def base(ax, title=None, scale=True, labels=True):
    ax.plot(basin[:,0],basin[:,1],color='k',lw=0.6)
    ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(aspect)
    ax.set_xticks([-77.8,-77.5,-77.2]); ax.set_yticks([-10.0,-9.5,-9.0])
    ax.set_xticklabels(['77.8°W','77.5°W','77.2°W'] if labels else []); ax.set_yticklabels(['10.0°S','9.5°S','9.0°S'] if labels else [])
    ax.tick_params(length=2,labelsize=6.5)
    if title: ax.set_title(title,fontsize=8.5)
    if scale:
        km=20; dlon=km/(111.32*np.cos(np.radians(lat0))); x0,y0=W+0.03,S+0.06
        ax.plot([x0,x0+dlon],[y0,y0],color='k',lw=2); ax.text(x0+dlon/2,y0+0.025,f'{km} km',ha='center',fontsize=6.5)
        ax.annotate('N',xy=(E-0.07,N-0.06),xytext=(E-0.07,N-0.17),ha='center',fontsize=8,fontweight='bold',arrowprops=dict(arrowstyle='-|>',color='k',lw=1))
CN=['0–7 Low','7–10 Slight','10–50 Moderate','50–100 Medium','100–150 Considerable','> 150 High']
CC=['#8c6d31','#b15ce6','#f4a582','#4daf4a','#ff9f1c','#d7191c']
# Fig 1 glaciers
fig,ax=plt.subplots(figsize=(3.6,5.6))
ax.imshow(img('red_drenaje'),extent=ext,interpolation='nearest',zorder=1,alpha=0.8)
ax.imshow(img('glaciar_actual'),extent=ext,interpolation='nearest',zorder=2)
ax.imshow(img('desglaciado'),extent=ext,interpolation='nearest',zorder=3)
base(ax)
ox,oy=m['outlet']; ax.plot(ox,oy,'r^',ms=6,zorder=5); ax.text(ox+0.04,oy,'La Balsa',fontsize=7,va='center')
ax.legend(handles=[Patch(color='#00bfff',label='Glacier 2023–2024 (Sentinel-2)'),Patch(color='#ff00ff',label='Deglaciated since 1987–1991 (Landsat 5)'),Line2D([],[],color='#1565c0',lw=1.5,label='Drainage network (≥ 30 km²)'),Line2D([],[],color='k',lw=0.8,label='Basin boundary')],loc='lower left',fontsize=6,frameon=True,bbox_to_anchor=(0.0,0.1))
fig.tight_layout(); fig.savefig(f'{OUT}/Figure1_study_area.png',dpi=300); plt.close(fig)
# Fig 2 factors
specs=[('factor_R','R factor (R₄)','MJ mm ha⁻¹ h⁻¹ yr⁻¹',0,3000,['#ffffcc','#41b6c4','#253494']),
       ('factor_K','K factor','t ha h ha⁻¹ MJ⁻¹ mm⁻¹',0.01,0.04,['#ffffd4','#fe9929','#993404']),
       ('factor_LS','LS factor (RUSLE)','dimensionless',0,20,['#f7f7f7','#636363','#000000']),
       ('factor_C','C factor (NDVI)','dimensionless',0,1,['#1a9850','#fee08b','#d73027'])]
fig,axs=plt.subplots(1,4,figsize=(7.2,3.6))
for i,(n,t,u,lo,hi,pal) in enumerate(specs):
    ax=axs[i]; ax.imshow(img(n),extent=ext,interpolation='nearest'); base(ax,t,scale=(i==0),labels=(i==0))
    cb=fig.colorbar(ScalarMappable(Normalize(lo,hi),LinearSegmentedColormap.from_list('p',pal)),ax=ax,orientation='horizontal',fraction=0.05,pad=0.06)
    cb.set_label(u,fontsize=6); cb.ax.tick_params(labelsize=5.5)
fig.tight_layout(); fig.savefig(f'{OUT}/Figure2_factors.png',dpi=300); plt.close(fig)
# Fig 3 models
fig,axs=plt.subplots(1,3,figsize=(7.2,4.2))
for i,(n,t) in enumerate([('erosion_USLE','USLE'),('erosion_RUSLE','RUSLE (reference)'),('erosion_RUSLE_MB','RUSLE with Moore–Burch LS')]):
    ax=axs[i]; ax.imshow(img(n),extent=ext,interpolation='nearest'); base(ax,t,scale=(i==0),labels=(i==0))
fig.legend(handles=[Patch(color=c,label=l) for c,l in zip(CC,CN)],loc='lower center',ncol=6,fontsize=6.5,title='Current potential soil loss (t ha⁻¹ yr⁻¹)',title_fontsize=7,frameon=False)
fig.tight_layout(rect=[0,0.08,1,1]); fig.savefig(f'{OUT}/Figure3_models.png',dpi=300); plt.close(fig)
# Fig 4 mining legacies
fig,ax=plt.subplots(figsize=(4.2,5.8))
ax.imshow(img('erosion_RUSLE'),extent=ext,interpolation='nearest',alpha=0.9)
P=m['pasivos']; xy=np.array([p['coords'][:2] for p in P]); h=np.array([p.get('hot') if p.get('hot') is not None else 0 for p in P])
cols=np.where(h>=0.3,'#000000',np.where(h>=0.1,'#4d4d4d','#ffffff'))
ax.scatter(xy[:,0],xy[:,1],s=9,c=cols,edgecolors='k',linewidths=0.4,zorder=5)
base(ax)
h1=[Patch(color=c,label=l) for c,l in zip(CC,CN)]
h2=[Line2D([],[],marker='o',ls='',mfc='#000000',mec='k',ms=4,label='≥ 30 % of 500 m buffer in hotspots'),Line2D([],[],marker='o',ls='',mfc='#4d4d4d',mec='k',ms=4,label='10–30 %'),Line2D([],[],marker='o',ls='',mfc='#ffffff',mec='k',ms=4,label='< 10 %')]
l1=ax.legend(handles=h1,loc='lower left',fontsize=5.8,title='RUSLE (t ha⁻¹ yr⁻¹)',title_fontsize=6.2,bbox_to_anchor=(0,0.1)); ax.add_artist(l1)
ax.legend(handles=h2,loc='upper right',bbox_to_anchor=(1.0,0.86),fontsize=5.8,title='Mining legacies (n = 399)',title_fontsize=6.2)
fig.tight_layout(); fig.savefig(f'{OUT}/Figure4_mining_legacies.png',dpi=300); plt.close(fig)
print('ok')
