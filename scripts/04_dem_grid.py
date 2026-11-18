"""Copernicus GLO-90 elevation and slope on the grid of the exported layers (results/dem_grid.npz).
Run from the repository root after 01_decode_layers.py. Tiles are downloaded from the public AWS bucket if missing."""
import os, urllib.request
import rasterio, numpy as np, json
os.makedirs('data/dem', exist_ok=True)
for t in ['S09_00_W078_00','S10_00_W078_00','S11_00_W078_00']:
    f=f'data/dem/cop_{t}.tif'
    if not os.path.exists(f):
        urllib.request.urlretrieve(f'https://copernicus-dem-90m.s3.amazonaws.com/Copernicus_DSM_COG_30_{t}_DEM/Copernicus_DSM_COG_30_{t}_DEM.tif', f)
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
m=json.load(open('data/layers/meta.json')); W,S,E,N=m['bounds']
srcs=[rasterio.open(f'data/dem/cop_{t}.tif') for t in ['S09_00_W078_00','S10_00_W078_00','S11_00_W078_00']]
mos,tr=merge(srcs)
g=np.load('results/grid.npz'); H,Wd=g['valid'].shape
dst=np.zeros((H,Wd),np.float32); dtr=from_bounds(W,S,E,N,Wd,H)
reproject(mos[0],dst,src_transform=tr,src_crs=srcs[0].crs,dst_transform=dtr,dst_crs=srcs[0].crs,resampling=Resampling.average)
# slope from native 90 m mosaic, then average to grid
z=mos[0].astype(np.float64); res=tr.a; lat=np.radians((N+S)/2)
dy=abs(tr.e)*111320; dx=res*111320*np.cos(lat)
gy,gx=np.gradient(z,dy,dx); sl=np.degrees(np.arctan(np.hypot(gx,gy))).astype(np.float32)
sdst=np.zeros((H,Wd),np.float32)
reproject(sl,sdst,src_transform=tr,src_crs=srcs[0].crs,dst_transform=dtr,dst_crs=srcs[0].crs,resampling=Resampling.average)
np.savez_compressed('results/dem_grid.npz',z=dst,slope=sdst)
v=g['valid']; print('elev valid mean',dst[v].mean(), 'min',dst[v].min(),'max',dst[v].max(),'slope mean',sdst[v].mean())
