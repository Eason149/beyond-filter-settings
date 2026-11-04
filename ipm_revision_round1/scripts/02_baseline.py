from common import *
import statsmodels.formula.api as smf
import warnings
b=behavior();rows=[];checks=[]
for fold,subjects in enumerate(folds(),1):
    te=b.Participant.isin(subjects)
    for model,extra in [('factor',[]),('factor_plus_A2',['A2'])]:
        pred=fit_predict(b[~te],b[te],extra)
        rows.append(dict(model=model,fold=fold,**scores(b.loc[te,'Rating'],pred)))
r=pd.DataFrame(rows);save(r,'baseline/participant_cv.csv')
s=r.groupby('model',sort=False)[['r2','rmse','mae']].mean().reset_index();save(s,'baseline/participant_cv_summary.csv')
old=pd.read_csv(REPO/'results/stage_2_7/a2_behavior_cv_summary.csv')
for _,row in s.iterrows():
    a=old[old.model.eq(row.model)].iloc[0]
    for k in ['r2','rmse','mae']:checks.append(dict(result=row.model,quantity=k,reproduced=row[k],archived=a[k],difference=row[k]-a[k],tolerance=1e-7))
print('Participant CV reproduced',flush=True)
with warnings.catch_warnings(record=True) as warns:
    fit=smf.mixedlm('Rating ~ RatingDimension_Beauty*(G_new+A2+I)+C(Identity)',b,groups=b.Participant,re_formula='~RatingDimension_Beauty').fit(reml=False,method='lbfgs',maxiter=1000)
coef=pd.DataFrame(dict(term=fit.params.index,estimate=fit.params.values,se=fit.bse.values,p=fit.pvalues.values))
save(coef,'baseline/joint_GAI_coefficients.csv')
(OUT/'baseline/joint_GAI_diagnostics.json').write_text(json.dumps(dict(converged=bool(fit.converged),aic=fit.aic,bic=fit.bic,llf=fit.llf,warnings=[str(w.message) for w in warns]),indent=2),encoding='utf-8')
old=pd.read_csv(REPO/'results/stage_2_6/behavior_frozen_metric_models.csv');old=old[old.model.eq('A2')]
for _,row in coef.iterrows():
    a=old[old.term.eq(row.term)].iloc[0]
    for k,tol in [('estimate',1e-7),('se',1e-4)]:checks.append(dict(result='joint_GAI:'+row.term,quantity=k,reproduced=row[k],archived=a[k],difference=row[k]-a[k],tolerance=tol))
print('Joint GAI mixed model reproduced',flush=True)
d=eeg()
for descriptor,seed,source,name in [(False,20261828,'results/stage_2_8/predefined_window_group_stats.csv','operation_eeg'),(True,20260904,'results/stage_2_9/gai_to_eeg_group_stats.csv','descriptor_eeg')]:
    bet,diag=eeg_betas(d,descriptor);g=group(bet,seed)
    save(bet,f'baseline/{name}_subject_betas.csv');save(g,f'baseline/{name}_group_stats.csv');save(diag,f'baseline/{name}_diagnostics.csv')
    old=pd.read_csv(REPO/source)
    if not descriptor:old=old.rename(columns={'factor':'metric'});old['window']=old.window.map(dict(zip(TAGS,W)))
    for _,row in g.iterrows():
        a=old[old.window.eq(row.window)&old.metric.eq(row.metric)].iloc[0]
        for k in ['mean_beta','se','p_maxT_familywise']:checks.append(dict(result=name+':'+row.window+':'+row.metric,quantity=k,reproduced=row[k],archived=a[k],difference=row[k]-a[k],tolerance=1e-7))
    print(name+' reproduced',flush=True)
checks=pd.DataFrame(checks);checks['pass']=checks.difference.abs().le(checks.tolerance);save(checks,'baseline/reproduction_comparison.csv')
print(checks.groupby('quantity')['pass'].agg(['sum','count']).to_string(),flush=True)
