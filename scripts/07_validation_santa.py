import json, numpy as np
from scipy import stats
a=json.load(open('results/annual_rusle_hy.json'))
meas={'2002-2003':3975971.8,'2003-2004':1364498.8,'2004-2005':3286383.7,'2005-2006':11112577.4,'2006-2007':5289617.6,'2007-2008':3280767.9,'2008-2009':8833249.8}
ks=list(meas); M=np.array([meas[k] for k in ks])/6815
Aha=np.array([a[k]['A'] for k in ks]); R=np.array([a[k]['R4'] for k in ks]); Pp=np.array([a[k]['P'] for k in ks]); MF=np.array([a[k]['MFI'] for k in ks])
area_ha=427750; B=4970.7
def y(sdr): return Aha*area_ha*sdr/B
res={}
for name,sdr in [('Vanoni',0.163),('Boyce',0.049),('USDA',0.222)]:
    Y=y(sdr); pb=100*(Y-M).sum()/M.sum(); nse=1-((Y-M)**2).sum()/((M-M.mean())**2).sum()
    res[name]=dict(Y=Y.round(0).tolist(),PBIAS=round(pb,1),NSE=round(nse,2))
r,p=stats.pearsonr(Aha,M); rs,ps=stats.spearmanr(Aha,M)
print('measured',M.round(0)); 
for k,v in res.items(): print(k,v)
print('pearson A vs meas r=%.2f p=%.3f ; spearman %.2f p=%.3f'%(r,p,rs,ps))
print('pearson P vs meas', stats.pearsonr(Pp,M)); print('pearson MFI vs meas',stats.pearsonr(MF,M))
print('CV measured %.0f%%  CV model %.0f%%'%(100*M.std(ddof=1)/M.mean(),100*Aha.std(ddof=1)/Aha.mean()))
print('ratio max/min meas %.1f model %.2f'%(M.max()/M.min(),Aha.max()/Aha.min()))
# SDR that matches mean
print('SDR matching 7yr mean:', round(M.mean()*B/(Aha.mean()*area_ha),3))
# power-law fit log(M)=a+b log(R)
b,a0,rr,pp,se=stats.linregress(np.log(R),np.log(M)); print('elasticity of measured yield to R4: b=%.1f (se %.1f) r=%.2f p=%.3f'%(b,se,rr,pp))
json.dump(dict(years=ks,measured_t_km2=M.round(1).tolist(),A_tha=Aha.round(1).tolist(),R4=R.round(0).tolist(),P=Pp.round(0).tolist(),MFI=MF.round(1).tolist(),sdr=res,pearson=[r,p],spearman=[rs,ps],elasticity=[b,se,rr,pp]),open('results/validation.json','w'),indent=1)
