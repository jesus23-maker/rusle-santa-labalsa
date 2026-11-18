"""Validation with annual RUSLE computed pixel by pixel in Earth Engine (30 m)."""
import json, numpy as np
from scipy import stats
ee = json.load(open('results/ee_results.json'))
meas_t = {'2002-2003':3975971.8,'2003-2004':1364498.8,'2004-2005':3286383.7,'2005-2006':11112577.4,'2006-2007':5289617.6,'2007-2008':3280767.9,'2008-2009':8833249.8}
ks = list(meas_t); M = np.array([meas_t[k] for k in ks]) / 6815
area_ha = 427750; B = 4970.7
out = {}
for tag in ['SepAug']:
    A = np.array([ee[f'{tag}_{k}']['A']['R4'] for k in ks]); P = np.array([ee[f'{tag}_{k}']['P']['P'] for k in ks])
    r = stats.pearsonr(A, M); rho = stats.spearmanr(A, M)
    res = dict(A=A.round(1).tolist(), P=P.round(0).tolist(), r=r[0], p=r[1], rho=rho[0], prho=rho[1], cv=100*A.std(ddof=1)/A.mean(), ratio=A.max()/A.min())
    for name, sdr in [('Vanoni', 0.163), ('Boyce', 0.049), ('USDA', 0.222)]:
        Y = A*area_ha*sdr/B
        res[name] = dict(Y=Y.round(0).tolist(), PBIAS=100*(Y-M).sum()/M.sum(), NSE=1-((Y-M)**2).sum()/((M-M.mean())**2).sum())
    res['sdr_match'] = M.mean()*B/(A.mean()*area_ha)
    pc = []
    for i in range(len(M)):
        k = np.arange(len(M)) != i; s = M[k].mean()*B/(A[k].mean()*area_ha); pc.append(A[i]*area_ha*s/B)
    pc = np.array(pc); res['loocv_const_NSE'] = 1-((pc-M)**2).sum()/((M-M.mean())**2).sum()
    Ycal = A*area_ha*res['sdr_match']/B; res['resid_sign'] = np.sign(M-Ycal).astype(int).tolist()
    o = np.argsort(A)[::-1][:2]; res['wet2_model_share'] = A[o].sum()/A.sum(); res['wet2_meas_share'] = M[o].sum()/M.sum(); res['wet2_years'] = [ks[j] for j in o]
    out[tag] = res
out['clim_A'] = ee['clim']['A']['R4']; out['clim_P'] = ee['clim']['P']['P']
out['measured'] = M.round(0).tolist(); out['years'] = ks
json.dump(out, open('results/ee_validation.json', 'w'), indent=1, default=float)
for tag in ['SepAug']:
    r = out[tag]
    print(tag, 'A', r['A'], 'r=%.2f p=%.3f rho=%.2f p=%.3f CV=%.1f%% ratio=%.2f' % (r['r'], r['p'], r['rho'], r['prho'], r['cv'], r['ratio']))
    print('   Vanoni PBIAS %.1f NSE %.2f | LOOCV const NSE %.2f | SDR match %.3f | resid %s | wet2 %s model %.2f meas %.2f' % (r['Vanoni']['PBIAS'], r['Vanoni']['NSE'], r['loocv_const_NSE'], r['sdr_match'], r['resid_sign'], r['wet2_years'], r['wet2_model_share'], r['wet2_meas_share']))
    print('   Y vanoni', r['Vanoni']['Y'])
print('clim', out['clim_A'])
