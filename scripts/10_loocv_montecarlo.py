import json, numpy as np
from scipy import stats
v=json.load(open('results/validation.json')); M=np.array(v['measured_t_km2']); R=np.array(v['R4']); A=np.array(v['A_tha'])
area_ha=427750; B=4970.7
# LOOCV of power-law yield model log M = a + b log R
pred=[]
for i in range(len(M)):
    k=np.arange(len(M))!=i; b,a,*_=stats.linregress(np.log(R[k]),np.log(M[k])); pred.append(np.exp(a+b*np.log(R[i])))
pred=np.array(pred); nse=1-((pred-M)**2).sum()/((M-M.mean())**2).sum(); pb=100*(pred-M).sum()/M.sum()
print('LOOCV erosivity-dependent model: NSE %.2f PBIAS %.1f'%(nse,pb), pred.round(0))
# constant calibrated SDR (LOOCV)
predc=[]
for i in range(len(M)):
    k=np.arange(len(M))!=i; sdr=M[k].mean()*B/(A[k].mean()*area_ha); predc.append(A[i]*area_ha*sdr/B)
predc=np.array(predc); print('LOOCV constant SDR: NSE %.2f'%(1-((predc-M)**2).sum()/((M-M.mean())**2).sum()))
# load concentration in two wettest years
o=np.argsort(R)[::-1][:2]; print('wettest 2 years', [v['years'][j] for j in o],'measured share %.2f model share %.2f'%(M[o].sum()/M.sum(), A[o].sum()/A.sum()))
# Monte Carlo
rng=np.random.default_rng(20261117); n=100000
fR=(1+rng.uniform(-0.15,0.11,n))**0.552        # CHIRPS bias range vs gauges, propagated through R4
fK=rng.uniform(1.0,1.65,n)                      # a=0.256 vs 0.0256
fC=rng.uniform(39.6/85.6,1.0,n)                 # tabulated vs NDVI C
Am=85.6*fR*fK*fC
lsdr=rng.uniform(np.log(0.049),np.log(0.222),n)  # log-uniform among Boyce..USDA
Y=Am*area_ha*np.exp(lsdr)/B
q=lambda x:np.percentile(x,[5,50,95]).round(1)
print('MC mean A t/ha/yr p5/50/95',q(Am)); print('MC yield p5/50/95',q(Y))
print('P(yield within measured 2002-09 range 200-1631) = %.2f ; P(within 779±162) = %.2f'%(((Y>=200)&(Y<=1631)).mean(),((Y>=617)&(Y<=941)).mean()))
# variance decomposition (first-order, log space)
L=np.log(Y); 
for nm,x in [('R',np.log(fR)),('K',np.log(fK)),('C',np.log(fC)),('SDR',lsdr)]: print(nm,'share of log-variance %.2f'%(np.var(x)/np.var(L)))
json.dump(dict(loocv_power=dict(NSE=nse,PBIAS=pb,pred=pred.tolist()),loocv_const_NSE=float(1-((predc-M)**2).sum()/((M-M.mean())**2).sum()),
 wet2_share_meas=float(M[o].sum()/M.sum()),wet2_share_model=float(A[o].sum()/A.sum()),
 mc_A=q(Am).tolist(),mc_Y=q(Y).tolist(),p_range=float(((Y>=200)&(Y<=1631)).mean()),p_ci=float(((Y>=617)&(Y<=941)).mean())),open('results/loocv_mc.json','w'),indent=1)
