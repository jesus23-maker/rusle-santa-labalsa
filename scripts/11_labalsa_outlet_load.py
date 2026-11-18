import numpy as np, json
Q={1975:[116.52,133.89,276.09,148.68,87.61,42.58,29.25,30.30,33.81,49.46,None,66.90,89.68],
1976:[114.91,157.77,171.88,94.89,47.35,33.79,23.59,22.11,25.84,44.62,None,63.98,70.89],
1977:[114.82,184.14,185.45,84.57,52.92,34.40,29.40,34.71,36.28,45.11,None,93.62,82.12],
1978:[96.77,136.92,109.67,81.17,52.41,36.13,31.96,31.47,41.21,42.54,76.48,89.44,68.86],
1979:[92.28,144.69,256.70,110.46,60.40,42.13,34.78,34.86,38.26,48.50,63.06,71.44,83.13],
1980:[90.09,86.71,82.16,68.27,40.29,35.83,28.44,31.23,46.57,63.88,78.08,116.54,64.01]}
for y,v in Q.items():
    if v[10] is None: v[10]=round(12*v[12]-sum(x for x in v[:12] if x is not None),2)
    print(y,'Nov',v[10],'check annual',round(np.mean(v[:12]),2),v[12])
M=np.array([Q[y][:12] for y in Q]).mean(0)
C=np.array([0.826,2.054,2.676,1.326]+[0.233]*8)   # g/L = kg/m3
days=np.array([31,28.25,31,30,31,30,31,31,30,31,30,31])
load=(C*M*days*86400/1000).sum()  # t/yr
print('mean monthly Q 1975-80',M.round(1),'annual',M.mean().round(1))
print('load Mt/yr %.2f'%(load/1e6), 'wet Jan-Apr share %.2f'%((C*M*days)[:4].sum()/(C*M*days).sum()))
for A in [4260,4768,4970.7]: print('yield over',A,'km2: %.0f t/km2/yr'%(load/A))
# flow-weighted concentration
print('flow-weighted C %.3f g/L'%((C*M*days).sum()/(M*days).sum()))
# per-year loads using monthly C climatology x each year's Q (variability from Q only)
yl={y:float((C*np.array(Q[y][:12])*days*86400/1000).sum()) for y in Q}
print({y:round(v/1e6,2) for y,v in yl.items()})
# sensitivity: max monthly C replacing mean is not used; use +-: dry-season C 0.12 (1997 week) instead of 0.233
C2=C.copy(); C2[4:]=0.12; print('alt dry C 0.12 -> Mt %.2f'%((C2*M*days*86400/1000).sum()/1e6))
json.dump(dict(Q_monthly_mean_1975_1980=M.round(2).tolist(),load_t_yr=load,yield_4260=load/4260,yield_4768=load/4768,yield_4970=load/4970.7,annual_loads=yl,Q_by_year=Q),open('results/labalsa_outlet_load.json','w'),indent=1)
