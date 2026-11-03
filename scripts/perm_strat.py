import json, numpy as np
exec(open('scripts/02_permutation_test.py').read().split('vr,vc=np.nonzero(valid)')[0])
d=np.load('results/dem_grid.npz'); z=d['z']; sl=d['slope']
eb=np.digitize(z,[3000,4000,4800]); sb=np.digitize(sl,[15,25,35])
strat=eb*4+sb
rng=np.random.default_rng(20261102)
vr,vc=np.nonzero(valid); vs=strat[vr,vc]
pools={s:np.where(vs==s)[0] for s in np.unique(vs)}
pst=np.array([strat[r,c] for r,c in [rc(*p['coords'][:2]) for p in keep]])
# precompute buffer stats for every valid pixel lazily via cache
cache={}
def b(i):
    if i not in cache: cache[i]=buf(vr[i],vc[i])
    return cache[i]
def rep_once():
    out=[]
    for s in pst:
        pool=pools.get(s)
        while True:
            i=pool[rng.integers(len(pool))]; x=b(i)
            if x: out.append(x); break
    return np.array(out).mean(0)
N3=np.array([rep_once() for _ in range(5000)])
res={}
for j,lab in [(1,'cls'),(0,'hot')]:
    o=obs[:,j].mean(); p=(np.sum(N3[:,j]>=o)+1)/(len(N3)+1)
    res[lab]=dict(obs=round(float(o),4),null=round(float(N3[:,j].mean()),4),lo=round(float(np.percentile(N3[:,j],2.5)),4),hi=round(float(np.percentile(N3[:,j],97.5)),4),p=round(float(p),4))
# describe strata of pasivos vs basin
import collections
pe=collections.Counter(eb[[rc(*p['coords'][:2])[0] for p in keep],[rc(*p['coords'][:2])[1] for p in keep]].tolist())
ps=sl[[rc(*p['coords'][:2])[0] for p in keep],[rc(*p['coords'][:2])[1] for p in keep]]
res['pasivos_slope_mean']=round(float(ps.mean()),2); res['basin_slope_mean']=round(float(sl[valid].mean()),2)
res['pasivos_elev_mean']=round(float(z[[rc(*p['coords'][:2])[0] for p in keep],[rc(*p['coords'][:2])[1] for p in keep]].mean()),0); res['basin_elev_mean']=round(float(z[valid].mean()),0)
res['pasivos_by_elevband']=dict(pe)
print(json.dumps(res,indent=1)); json.dump(res,open('results/perm_strat_results.json','w'),indent=1)
