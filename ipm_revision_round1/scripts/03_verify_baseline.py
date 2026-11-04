from common import *
c=pd.read_csv(OUT/'baseline/reproduction_comparison.csv')
old=pd.read_csv(REPO/'results/chb_stimulus_robustness/stimulus_aware_A2_increment.csv').iloc[0]
new=pd.read_csv(OUT/'baseline/crossed/stimulus_aware_A2_increment.csv').iloc[0]
rows=[]
for k in ['LR_chi2','delta_AIC','delta_BIC']:
    rows.append(dict(result='crossed_A2_increment',quantity=k,reproduced=new[k],archived=old[k],difference=new[k]-old[k],tolerance=.01))
a=pd.read_csv(REPO/'results/chb_stimulus_robustness/stimulus_aware_coefficients.csv')
b=pd.read_csv(OUT/'baseline/crossed/stimulus_aware_coefficients.csv')
for _,r in b.iterrows():
    o=a[a.model.eq(r.model)&a.term.eq(r.term)].iloc[0]
    for k,tol in [('estimate',1e-7),('se',1e-4)]:rows.append(dict(result=r.model+':'+r.term,quantity=k,reproduced=r[k],archived=o[k],difference=r[k]-o[k],tolerance=tol))
r=pd.DataFrame(rows);r['pass']=r.difference.abs().le(r.tolerance);c=pd.concat([c,r],ignore_index=True);save(c,'baseline/reproduction_comparison.csv')
summary=pd.read_csv(OUT/'baseline/participant_cv_summary.csv')
report=f'''# Baseline reproduction

Executed from the archived behavioral and EEG feature inputs, using the fixed model specifications. Source EEG acquisition/preprocessing and MediaPipe/embedding extraction are audited from archived scripts; they are not rerun in this round. A2 algebra was independently verified from its stored components.

Checks within their prespecified tolerances: {c['pass'].sum()} / {len(c)}. Full comparison: `baseline/reproduction_comparison.csv`. Coefficients/scores tolerance 1e-7, SE 1e-4, optimizer-dependent likelihood differences 0.01. All input counts are asserted before estimation.

Participant-held-out performance:
```
{summary.to_string(index=False)}
```

Archived stimulus-intercept crossed model: reproduced LR={new.LR_chi2:.9f}, delta AIC={new.delta_AIC:.9f}, delta BIC={new.delta_BIC:.9f}. Fits and convergence diagnostics: `baseline/crossed/`.

Predefined operation and descriptor EEG coefficients and max-|t| p-values are in `baseline/*eeg_group_stats.csv`; families, permutation count, seed and subject order were preserved. No window or channel searches were executed.

Discrepancies:
```
{c[~c['pass']].to_string(index=False) if (~c['pass']).any() else 'None outside tolerance.'}
```

This baseline retains the original 15 zero-valued rating entries for exact comparison; their exclusion is tested separately. New analyses may now proceed under `analysis_plan.md`. Reproduction establishes numerical lineage, while the additional stimulus-dimension-slope sensitivity examines inference beyond the original covariance specification.
'''
(OUT/'baseline_reproduction.md').write_text(report,encoding='utf-8');print(report)
assert c['pass'].all(), 'Investigate discrepancy before new models'
