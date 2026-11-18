import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
v=json.load(open('results/validation.json')); M=np.array(v['measured_t_km2']); Y=np.array(v['sdr']['Vanoni']['Y'])
lo=np.array(v['sdr']['Boyce']['Y']); hi=np.array(v['sdr']['USDA']['Y']); R=np.array(v['R4']); yrs=[s.replace('-20','–') for s in v['years']]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
fig,ax=plt.subplots(1,3,figsize=(11,3.4),gridspec_kw={'width_ratios':[1.5,1,1]})
x=np.arange(len(M)); w=0.38
ax[0].bar(x-w/2,M,w,color='#2b6c8f',label='Measured, Santa station (Morera, 2010)')
ax[0].errorbar(x-w/2,M,yerr=[M*0.38,M*0.85],fmt='none',ecolor='#14384b',lw=0.8,capsize=2)
ax[0].bar(x+w/2,Y,w,color='#e08a3c',label='RUSLE × SDR (Vanoni)')
ax[0].errorbar(x+w/2,Y,yerr=[Y-lo,hi-Y],fmt='none',ecolor='#7a4515',lw=0.8,capsize=2)
ax[0].set_xticks(x); ax[0].set_xticklabels(yrs,rotation=35,ha='right'); ax[0].set_ylabel('Suspended sediment yield (t km$^{-2}$ yr$^{-1}$)')
ax[0].legend(frameon=False,fontsize=7.5,loc='upper left'); ax[0].set_title('(a) Hydrological years 2002–2009',loc='left',fontsize=9)
ax[1].scatter(Y,M,c='#2b6c8f',s=28,zorder=3)
for i,t in enumerate(yrs): ax[1].annotate(t[2:4]+'–'+t[-2:],(Y[i],M[i]),fontsize=6.5,xytext=(3,2),textcoords='offset points')
l=[0,2000]; ax[1].plot(l,l,'--',c='grey',lw=0.8); ax[1].set_xlim(900,1500); ax[1].set_ylim(0,1800)
ax[1].set_xlabel('Modelled (t km$^{-2}$ yr$^{-1}$)'); ax[1].set_ylabel('Measured'); ax[1].set_title('(b) r = %.2f, ρ = %.2f'%(v['pearson'][0],v['spearman'][0]),loc='left',fontsize=9)
ax[2].scatter(R,M,c='#2b6c8f',s=28,zorder=3); b,se,rr,pp=v['elasticity']
rg=np.linspace(R.min()*0.97,R.max()*1.03,50); a=np.exp(np.mean(np.log(M))-b*np.mean(np.log(R))); ax[2].plot(rg,a*rg**b,c='#2b6c8f',lw=1)
ax[2].plot(rg,M.mean()*(rg/R.mean()),c='#e08a3c',lw=1,ls='--')
ax[2].set_xlabel('Basin-mean R$_4$ (MJ mm ha$^{-1}$ h$^{-1}$ yr$^{-1}$)'); ax[2].set_ylabel('Measured yield')
ax[2].set_title('(c) Yield ∝ R$^{%.1f}$ (RUSLE: R$^{1}$)'%b,loc='left',fontsize=9)
for a_ in ax: a_.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('results/Figure5_validation.png',dpi=300)
