import numpy as np, json, rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.features import rasterize
from PIL import Image
from scipy import ndimage
code=''
m=json.load(open(code+'data/layers/meta.json')); W,S,E,N=m['bounds']
srcs=[rasterio.open(f'data/dem/cop_{t}.tif') for t in ['S09_00_W078_00','S10_00_W078_00','S11_00_W078_00']]
mos,tr=merge(srcs,bounds=(W-0.02,S-0.02,E+0.02,N+0.02))
crs='EPSG:32718'; res=90.0
dtr,dw,dh=calculate_default_transform('EPSG:4326',crs,mos.shape[2],mos.shape[1],W-0.02,S-0.02,E+0.02,N+0.02,resolution=res)
z=np.zeros((dh,dw),np.float32); reproject(mos[0],z,src_transform=tr,src_crs='EPSG:4326',dst_transform=dtr,dst_crs=crs,resampling=Resampling.bilinear)
# basin mask
from pyproj import Transformer
T=Transformer.from_crs('EPSG:4326',crs,always_xy=True)
poly=[T.transform(x,y) for x,y in m['basin']['coordinates'][0]]
basin=rasterize([({'type':'Polygon','coordinates':[poly]},1)],out_shape=z.shape,transform=dtr).astype(bool)
# erosion class grid -> A (t/ha/yr) in UTM
img=np.array(Image.open(code+'data/layers/erosion_RUSLE.png').convert('RGBA'))
COLS=np.array([[0x8c,0x6d,0x31],[0xb1,0x5c,0xe6],[0xf4,0xa5,0x82],[0x4d,0xaf,0x4a],[0xff,0x9f,0x1c],[0xd7,0x19,0x1c]])
valid=img[...,3]>0; cls=((img[...,None,:3].astype(int)-COLS[None,None])**2).sum(-1).argmin(-1)
mids=np.array([3.5,8.5,30,75,125,225.0]); Ag=np.where(valid,mids[cls],np.nan).astype(np.float32)
from rasterio.transform import from_bounds
gtr=from_bounds(W,S,E,N,Ag.shape[1],Ag.shape[0])
A=np.full(z.shape,np.nan,np.float32); reproject(Ag,A,src_transform=gtr,src_crs='EPSG:4326',dst_transform=dtr,dst_crs=crs,resampling=Resampling.nearest,src_nodata=np.nan,dst_nodata=np.nan)
A[~basin]=np.nan
# flow routing with pysheds (D8)
np.in1d=np.isin
from pysheds.sview import Raster, ViewFinder
from pysheds.grid import Grid
import pyproj; vf=ViewFinder(affine=dtr,shape=z.shape,crs=pyproj.Proj('EPSG:32718'),nodata=-9999.)
dem=Raster(np.where(basin|True,z,-9999).astype('f8'),viewfinder=vf)
g=Grid(viewfinder=vf)
d=g.fill_pits(dem); d=g.fill_depressions(d); d=g.resolve_flats(d)
dirmap=(64,128,1,2,4,8,16,32)
fdir=np.asarray(g.flowdir(d,dirmap=dirmap)); acc=np.asarray(g.accumulation(g.flowdir(d,dirmap=dirmap),dirmap=dirmap)).astype('f8')
H,Wd=z.shape
# downstream index
off={64:(-1,0),128:(-1,1),1:(0,1),2:(1,1),4:(1,0),8:(1,-1),16:(0,-1),32:(-1,-1)}
rr,cc=np.indices(z.shape); dr=np.zeros_like(rr); dc=np.zeros_like(cc)
for k,(a,b) in off.items(): msk=fdir==k; dr[msk]=a; dc[msk]=b
nr=rr+dr; nc=cc+dc; hasd=(fdir>0)&(nr>=0)&(nr<H)&(nc>=0)&(nc<Wd)
down=np.where(hasd,nr*Wd+nc,-1).ravel()
dist=np.where((dr!=0)&(dc!=0),res*np.sqrt(2),res).ravel()
# slope and weight (Cavalli 2013 roughness)
gy,gx=np.gradient(d.astype('f8'),res); Sl=np.clip(np.hypot(gx,gy),0.005,1.0)
resid=z-ndimage.uniform_filter(z,5); RI=np.sqrt(np.maximum(ndimage.uniform_filter(resid**2,5)-ndimage.uniform_filter(resid,5)**2,0))
Wt=np.clip(1-RI/RI[basin].max(),0.001,1)
# upslope sums in increasing acc order
order=np.argsort(acc.ravel(),kind='stable')
sw=Wt.ravel().copy(); ss=Sl.ravel().copy(); n=np.ones(H*Wd)
dn=down
for i in order:
    j=dn[i]
    if j>=0: sw[j]+=sw[i]; ss[j]+=ss[i]; n[j]+=n[i]
Dup=(sw/n)*(ss/n)*np.sqrt(n*res*res)
# targets: stream cells (contributing area >= 5 km2)
stream=(acc.ravel()*res*res>=5e6)
cost=dist/(Wt.ravel()*Sl.ravel())
Ddn=np.zeros(H*Wd)
for i in order[::-1]:
    if stream[i]: Ddn[i]=0; continue
    j=dn[i]; Ddn[i]=cost[i]+(Ddn[j] if j>=0 else 0)
IC=np.log10(Dup/np.maximum(Ddn,res)).reshape(z.shape)
np.savez_compressed('results/ic_grid.npz',IC=IC.astype('f4'),A=A,basin=basin,stream=stream.reshape(z.shape),z=z)
ok=np.isfinite(A)&~stream.reshape(z.shape)
print('cells',ok.sum(),'IC mean %.2f p10 %.2f p90 %.2f'%(IC[ok].mean(),*np.percentile(IC[ok],[10,90])))
