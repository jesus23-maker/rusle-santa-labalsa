import json, numpy as np
rng=np.random.default_rng(20261030)
m=json.load(open('data/layers/meta.json')); W,S,E,N=m['bounds']
g=np.load('results/grid.npz'); cls=g['cls']; valid=g['valid']; hot=g['hot']
H,Wd=valid.shape; dlat=(N-S)/H; dlon=(E-W)/Wd
lat0=(N+S)/2; ky=111.32/ (dlat*111.32)  # px per degree lat (not used)
py_km=dlat*111.32; px_km=dlon*111.32*np.cos(np.radians(lat0))
R=0.5
ry=int(np.ceil(R/py_km))+1; rx=int(np.ceil(R/px_km))+1
oy,ox=np.mgrid[-ry:ry+1,-rx:rx+1]; disk=((oy*py_km)**2+(ox*px_km)**2)<=R**2
oy,ox=oy[disk],ox[disk]
def buf(r,c):
    rr=r+oy; cc=c+ox; k=(rr>=0)&(rr<H)&(cc>=0)&(cc<Wd); rr,cc=rr[k],cc[k]
    v=valid[rr,cc]
    if v.sum()==0: return None
    return hot[rr,cc][v].mean(), cls[rr,cc][v].mean()
def rc(lon,lat): return int((N-lat)/dlat), int((lon-W)/dlon)
P=[p for p in m['pasivos']]
obs=[]; keep=[]
for p in P:
    r,c=rc(*p['coords'][:2]); b=buf(r,c)
    if b: obs.append(b); keep.append(p)
obs=np.array(obs); n=len(obs)
units=np.array([str(p.get('exunidad')) for p in keep])
print('n',n,'obs hot',obs[:,0].mean().round(4),'obs cls',obs[:,1].mean().round(3),'EE hot',np.mean([p['hot'] for p in keep if p.get('hot') is not None]).round(4))
vr,vc=np.nonzero(valid)
# Null 1: random over basin
def rand_basin(k):
    out=[]
    while len(out)<k:
        i=rng.integers(len(vr)); b=buf(vr[i],vc[i])
        if b: out.append(b)
    return np.array(out)
N1=np.array([rand_basin(n).mean(0) for _ in range(5000)])
# Null 2: cluster-preserving displacement (whole ex-unidad shifted 2-10 km)
cent={u:None for u in sorted(set(units))}
prc=np.array([rc(*p['coords'][:2]) for p in keep])
def disp():
    res=[]
    for u in cent:
        idx=np.where(units==u)[0]
        for _ in range(200):
            dist=rng.uniform(2,10); ang=rng.uniform(0,2*np.pi)
            dr=int(round(-dist*np.sin(ang)/py_km)); dc=int(round(dist*np.cos(ang)/px_km))
            vals=[]; ok=True
            for i in idx:
                r,c=prc[i][0]+dr,prc[i][1]+dc
                if not(0<=r<H and 0<=c<Wd and valid[r,c]): ok=False;break
                b=buf(r,c)
                if b is None: ok=False;break
                vals.append(b)
            if ok: res+=vals; break
        else:
            res+= [tuple(obs[i]) for i in idx]
    return np.array(res).mean(0)
N2=np.array([disp() for _ in range(5000)])
def rep(name,Nd,j,lab):
    o=obs[:,j].mean(); mu=Nd[:,j].mean(); p=(np.sum(Nd[:,j]>=o)+1)/(len(Nd)+1)
    print(f'{name} {lab}: obs {o:.4f} null mean {mu:.4f} [{np.percentile(Nd[:,j],2.5):.4f},{np.percentile(Nd[:,j],97.5):.4f}] ratio {o/mu:.3f} p={p:.4f}')
    return dict(obs=o,null=mu,lo=np.percentile(Nd[:,j],2.5),hi=np.percentile(Nd[:,j],97.5),p=p)
out={}
for nm,Nd in [('cuenca',N1),('desplazamiento',N2)]:
    out[nm]={'hot':rep(nm,Nd,0,'hot'),'cls':rep(nm,Nd,1,'cls')}
# per-province with local null (point-level displacement)
prov=np.array([str(p.get('provincia')) for p in keep])
out['n']=n; out['units']=len(set(units))
json.dump(out,open('results/perm_results.json','w'),default=float,indent=1)
