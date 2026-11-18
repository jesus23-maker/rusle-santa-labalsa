import numpy as np, json
from PIL import Image
from matplotlib.path import Path
code=''
m=json.load(open(code+'data/layers/meta.json'))
W,S,E,N=m['bounds']
d=np.load('data/external/chirps_monthly_1995_2024.npz'); P=d['P']; yrs=d['years']; mos=d['months']; t=d['transform']
ny,nx=P.shape[1:]
# fine class grid
img=np.array(Image.open(code+'data/layers/erosion_RUSLE.png').convert('RGBA'))
COLS=np.array([[0x8c,0x6d,0x31],[0xb1,0x5c,0xe6],[0xf4,0xa5,0x82],[0x4d,0xaf,0x4a],[0xff,0x9f,0x1c],[0xd7,0x19,0x1c]])
valid=img[...,3]>0
cls=((img[...,None,:3].astype(int)-COLS[None,None])**2).sum(-1).argmin(-1)
pct=np.array([10.2,3.3,29.4,24.5,15.8,16.9])/100
mids=np.array([3.5,8.5,30,75,125,np.nan])
x=(85.6-np.nansum(pct[:5]*mids[:5]))/pct[5]; mids[5]=x
print('class repr >150:',round(x,1),' class shares in png:',np.round(np.bincount(cls[valid],minlength=6)/valid.sum(),3))
A=np.where(valid,mids[cls],0.0)
H,Wd=A.shape
lat=N-(np.arange(H)+.5)*(N-S)/H; lon=W+(np.arange(Wd)+.5)*(E-W)/Wd
LON,LAT=np.meshgrid(lon,lat)
ci=((LON-t[2])/t[0]).astype(int); ri=((LAT-t[5])/t[4]).astype(int)
wgt=np.zeros((ny,nx)); np.add.at(wgt,(ri[valid],ci[valid]),A[valid])
cnt=np.zeros((ny,nx)); np.add.at(cnt,(ri[valid],ci[valid]),1)
# basin mask at chirps res (any modelled pixel)
basin_cells=cnt>0
def R4(Pm):  # Pm: (12,ny,nx) monthly totals
    Pa=Pm.sum(0); mfi=(Pm**2).sum(0)/Pa; return 161*mfi**0.552, mfi, Pa
# climatology from mean monthly 1995-2024
clim=np.stack([P[mos==k].mean(0) for k in range(1,13)])
R4c,MFIc,Pc=R4(clim)
aw=lambda X:(X*cnt)[basin_cells].sum()/cnt[basin_cells].sum()
print('Clim: P %.0f  MFI %.1f  R4 %.0f (paper 778/96.8/1961)'%(aw(Pc),aw(MFIc),aw(R4c)))
def year_stats(sel):
    Pm=P[sel]; r,mf,pa=R4(Pm)
    ratio=(wgt*r/R4c)[basin_cells].sum()/wgt[basin_cells].sum()
    return dict(P=aw(pa),MFI=aw(mf),R4=aw(r),A=85.6*ratio)
out={}
# calendar 2017 check
s=year_stats(yrs==2017); print('2017 calendar: P %.0f MFI %.1f R4 %.0f A %.1f (paper 956/143.2/2433/106.1)'%(s['P'],s['MFI'],s['R4'],s['A']))
for y in range(1995,2024):
    sel=((yrs==y)&(mos>=9))|((yrs==y+1)&(mos<=8))
    out[f'{y}-{y+1}']=year_stats(sel)
json.dump(out,open('results/annual_rusle_hy.json','w'),indent=1)
for k in ['2002-2003','2003-2004','2004-2005','2005-2006','2006-2007','2007-2008','2008-2009']:
    v=out[k]; print(k,{a:round(b,1) for a,b in v.items()})
