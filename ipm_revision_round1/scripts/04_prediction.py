"""Frozen ablation, conditional paired uncertainty, identity transfer, influence."""
from common import *
import statsmodels.formula.api as smf
b=behavior();im=metrics();rows=[];pred=b[['Participant','Stimulus','Identity','Rating','RatingDimension_Beauty']].copy()
foldmap={s:k for k,f in enumerate(folds(),1) for s in f};pred['fold']=b.Participant.map(foldmap)
def scale(train,test,refit_a2=False):
    tr=train.copy();te=test.copy();u=tr.drop_duplicates('Stimulus')
    if refit_a2:
        rms=np.sqrt((u[COMP]**2).mean())
        z=(u[COMP]/rms).mean(axis=1);mu=z.mean();sd=z.std(ddof=0)
        for q in [tr,te]:q['A2']=((q[COMP]/rms).mean(axis=1)-mu)/sd
    for c in COMP:
        mu=u[c].mean();sd=u[c].std(ddof=0)
        for q in [tr,te]:q[c]=(q[c]-mu)/sd
    return tr,te
def result_rows(y,p,dim,baseinfo):
    return [dict(**baseinfo,dimension=label,n=int(mask.sum()),**scores(y[mask],p[mask])) for label,mask in [('pooled',np.ones(len(y),bool)),('unnaturalness',dim==0),('beauty',dim==1)]]
for k,subjects in enumerate(folds(),1):
    mask=b.Participant.isin(subjects);tr,te=scale(b[~mask],b[mask])
    for m,extra in MODELS.items():
        p=fit_predict(tr,te,extra);pred.loc[mask,m]=p
        rows+=result_rows(te.Rating.to_numpy(),p,te.RatingDimension_Beauty.to_numpy(),dict(model=m,fold=str(k),validation='participant'))
r=pd.DataFrame(rows)
mean=r.groupby(['model','dimension','validation'],sort=False)[['r2','rmse','mae']].mean().reset_index();mean['fold']='mean';mean['n']=np.where(mean.dimension.eq('pooled'),len(b),len(b)//2)
r=pd.concat([r,mean],ignore_index=True)
base=r[r.model.eq('M0')].set_index(['fold','dimension','validation'])
for m in ['r2','rmse','mae']:
    r['delta_'+m+'_vs_M0']=[row[m]-base.loc[(row.fold,row.dimension,row.validation),m] for _,row in r.iterrows()]
r['mse_reduction_pct_vs_M0']=[100*(1-row.rmse**2/base.loc[(row.fold,row.dimension,row.validation),'rmse']**2) for _,row in r.iterrows()]
for idx,row in r[r.fold.eq('mean')].iterrows():
    q=r[r.model.eq(row.model)&r.dimension.eq(row.dimension)&~r.fold.eq('mean')]
    q0=r[r.model.eq('M0')&r.dimension.eq(row.dimension)&~r.fold.eq('mean')]
    r.loc[idx,'mse_reduction_pct_vs_M0']=100*(1-np.mean(q.rmse**2)/np.mean(q0.rmse**2))
save(r,'a2_ablation_results.csv');save(pred,'tables/participant_oof_predictions.csv')
print('Ablation folds complete',flush=True)

# Evaluation-sample bootstrap: preserve paired predictions and fold estimand.
# Aggregate repeated outcomes exactly using count, sums and losses, never treating trials as clusters.
keys=['Participant','Stimulus','Identity','fold'];z=pred.copy();z['y2']=z.Rating**2
for m in MODELS:z[m+'_sse']=(z.Rating-z[m])**2;z[m+'_sae']=abs(z.Rating-z[m])
agg=z.groupby(keys,sort=True).agg(n=('Rating','size'),sy=('Rating','sum'),sy2=('y2','sum'),**{m+'_sse':(m+'_sse','sum') for m in MODELS},**{m+'_sae':(m+'_sae','sum') for m in MODELS}).reset_index()
subs=sorted(agg.Participant.unique(),key=lambda s:int(s[1:]));stims=sorted(agg.Stimulus.unique());si=agg.Participant.map({s:i for i,s in enumerate(subs)}).to_numpy();ii=agg.Stimulus.map({s:i for i,s in enumerate(stims)}).to_numpy()
idm={id:[stims.index(s) for s in agg.loc[agg.Identity.eq(id),'Stimulus'].unique()] for id in IDS}
statcols=['n','sy','sy2']+[m+'_sse' for m in MODELS]+[m+'_sae' for m in MODELS];a=agg[statcols].to_numpy();rng=np.random.default_rng(20261016);bs=[]
for unit in ['participant','stimulus_within_identity','participant_x_stimulus']:
    draws=[]
    for batch in range(10):
        ns=rng.multinomial(30,np.ones(30)/30,size=200) if 'participant' in unit else np.ones((200,30))
        ni=np.ones((200,64))
        if 'stimulus' in unit:
            for id,ix in idm.items():ni[:,ix]=rng.multinomial(len(ix),np.ones(len(ix))/len(ix),size=200)
        weights=ns[:,si]*ni[:,ii];fs=[]
        for k in range(1,6):
            select=agg.fold.eq(k).to_numpy();tot=weights[:,select]@a[select];n=tot[:,0,None];sse=tot[:,3:9];sae=tot[:,9:15];sst=tot[:,2,None]-tot[:,1,None]**2/n
            fs.append(np.stack([1-sse/sst,np.sqrt(sse/n),sae/n],axis=2))
        draws.append(np.mean(fs,axis=0))
    draws=np.concatenate(draws)
    for lhs,rhs in [(m,'M0') for m in list(MODELS)[1:]]+[('M4',m) for m in ['M1','M2','M3','M5']]:
        delta=draws[:,list(MODELS).index(lhs),:]-draws[:,list(MODELS).index(rhs),:]
        l=r[(r.fold=='mean')&(r.dimension=='pooled')&(r.model==lhs)].iloc[0];rr=r[(r.fold=='mean')&(r.dimension=='pooled')&(r.model==rhs)].iloc[0]
        for j,metric in enumerate(['r2','rmse','mae']):
            finite=delta[np.isfinite(delta[:,j]),j];lo,hi=np.quantile(finite,[.025,.975])
            bs.append(dict(bootstrap_unit=unit,comparison=lhs+' minus '+rhs,metric=metric,estimate=l[metric]-rr[metric],ci95_low=lo,ci95_high=hi,n_bootstrap=len(finite),n_attempted=2000,n_discarded_empty_fold=2000-len(finite),seed=20261016,conditional_on_fitted_models=True))
save(pd.DataFrame(bs),'tables/ablation_paired_uncertainty.csv')
print('Paired participant/image/product bootstrap complete',flush=True)

rows=[];pr=[]
for id in IDS:
    tr,te=scale(b[~b.Identity.eq(id)],b[b.Identity.eq(id)],refit_a2=True)
    q=te[['Participant','Stimulus','Identity','Rating','RatingDimension_Beauty']].copy()
    for m,extra in MODELS.items():
        p=fit_predict(tr,te,extra,identity=False);q[m]=p
        rows+=result_rows(te.Rating.to_numpy(),p,te.RatingDimension_Beauty.to_numpy(),dict(model=m,held_out_identity=id,n_train=len(tr)))
    pr.append(q)
r=pd.DataFrame(rows);base=r[r.model.eq('M0')].set_index(['held_out_identity','dimension'])
for m in ['r2','rmse','mae']:r['delta_'+m+'_vs_M0']=[row[m]-base.loc[(row.held_out_identity,row.dimension),m] for _,row in r.iterrows()]
save(r,'identity_holdout_results.csv');save(pd.concat(pr),'tables/identity_oof_predictions.csv')
print('Identity holdout complete',flush=True)

# Zero-entry sensitivity uses the same participant memberships and no imputation.
clean=b[b.Rating.between(1,7)].copy();zr=[]
for k,subjects in enumerate(folds(),1):
    tr,te=scale(clean[~clean.Participant.isin(subjects)],clean[clean.Participant.isin(subjects)])
    for m in ['M0','M4']:
        p=fit_predict(tr,te,MODELS[m]);zr.append(dict(model=m,fold=k,n_train=len(tr),n_test=len(te),**scores(te.Rating,p)))
save(pd.DataFrame(zr),'tables/zero_exclusion_cv.csv')
# Coefficient influence is a descriptive fixed-effect refit, not independent tests.
ir=[]
for unit,col in [('participant','Participant'),('stimulus','Stimulus'),('identity','Identity')]:
    for omit in sorted(b[col].unique()):
        q=b[~b[col].eq(omit)];fit=smf.ols(formula(['A2']),q).fit()
        v=fit.params['A2'];inter=fit.params['RatingDimension_Beauty:A2']
        ir.append(dict(omission_unit=unit,omitted=omit,n=len(q),unnaturalness_slope=v,beauty_slope=v+inter,dimension_interaction=inter))
save(pd.DataFrame(ir),'tables/behavior_coefficient_influence.csv')
save(im.loc[im[M].abs().gt(3).any(axis=1),['identity','picture_filename']+M],'tables/extreme_descriptor_images.csv')
print('Zero sensitivity and all participant/stimulus/identity coefficient influences complete',flush=True)
