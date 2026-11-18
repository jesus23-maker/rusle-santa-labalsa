import numpy as np, json
from scipy.optimize import brentq
g=np.load('results/ic_grid.npz'); IC=g['IC']; A=g['A']; st=g['stream']
ok=np.isfinite(A)&np.isfinite(IC)
a=A[ok]; ic=IC[ok]; s=st[ok]
def sdr(ic0,kb=2.0,smax=0.8): return np.where(s,smax,smax/(1+np.exp((ic0-ic)/kb)))
def eff(ic0,kb=2.0): x=sdr(ic0,kb); return (a*x).sum()/a.sum()
print('stream fraction of modelled cells %.3f'%s.mean())
print('InVEST defaults (IC0=0.5,kb=2): basin SDR = %.3f'%eff(0.5))
Aclim=85.6; area_ha=427750; B=4970.7
def yld(x,A_=Aclim): return A_*area_ha*x/B
print(' -> yield %.0f t/km2/yr (measured mean 779)'%yld(eff(0.5)))
for kb in [1.0,2.0,3.0]:
    ic0=brentq(lambda q: eff(q,kb)-0.104,-10,10); print('kb=%.0f calibrated IC0=%.2f'%(kb,ic0))
# spatial: fraction of delivered sediment from top-10% A cells vs fraction of eroded
x=sdr(0.5); p90=np.percentile(a,90)
print('hotspot share of gross erosion %.2f ; of delivered %.2f'%(a[a>=p90].sum()/a.sum(), (a*x)[a>=p90].sum()/(a*x).sum()))
json.dump(dict(SDR_invest_default=float(eff(0.5)),yield_invest_default=float(yld(eff(0.5))),IC_mean=float(ic.mean())),open('results/sdr_ic.json','w'),indent=1)
