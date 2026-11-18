import rasterio, numpy as np, json, os, sys
from rasterio.windows import from_bounds
from concurrent.futures import ThreadPoolExecutor
os.environ['GDAL_DISABLE_READDIR_ON_OPEN']='EMPTY_DIR'

W,S,E,N=-78.0,-10.25,-77.15,-8.75
U='/vsicurl/https://data.chc.ucsb.edu/products/CHIRPS/v2.0/global_monthly/cogs/chirps-v2.0.{y}.{m:02d}.cog'
def get(ym):
    y,m=ym
    for t in range(4):
        try:
            with rasterio.open(U.format(y=y,m=m)) as ds:
                w=from_bounds(W,S,E,N,ds.transform); a=ds.read(1,window=w).astype('f4'); tr=ds.window_transform(w)
                a[a<0]=np.nan; return ym,a,tuple(tr)[:6]
        except Exception as e: err=e
    raise err
yms=[(y,m) for y in range(1995,2025) for m in range(1,13)]
with ThreadPoolExecutor(12) as ex: res=list(ex.map(get,yms))
arr=np.stack([r[1] for r in res]); print(arr.shape, res[0][2])
np.savez_compressed('data/external/chirps_monthly_1995_2024.npz',P=arr,years=np.array([r[0][0] for r in res]),months=np.array([r[0][1] for r in res]),transform=np.array(res[0][2]))
