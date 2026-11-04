"""Read-only source loading and frozen analysis definitions for Round 1."""
from pathlib import Path
import os, hashlib, json, sys
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy import stats
import patsy
OUT=Path(__file__).resolve().parents[1]
REPO=OUT.parent
DATA=Path(os.environ.get('IPM_DATA_ROOT',REPO/'data')).resolve()
F=['FSlim','Eye','Mouth','Skin']
COMP=['a2_luminance','a2_chroma_js','a2_texture']
M=['G_new','A2','I']
IDS=['F_1','F_2','M_1','M_2']
SUB=[f's{i}' for i in range(1,31)]
SUB_EEG=[s for s in SUB if s not in ['s5','s18']]
W=['350-600 ms','600-1000 ms']
TAGS=['MiddleLate_350_600','Late_600_1000']
EC=[f'Centroparietal_{t}_Mean' for t in TAGS]
MODELS={'M0':[], 'M1':[COMP[0]],'M2':[COMP[1]],'M3':[COMP[2]],'M4':['A2'],'M5':COMP}
for p in ['tables','figures','logs','baseline','models']:(OUT/p).mkdir(exist_ok=True,parents=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(df,path):df.to_csv(OUT/path,index=False,encoding='utf-8-sig')
def metrics():
    if os.environ.get('IPM_ROUND1_REPLAY')=='1' or not (REPO/'results/stage_2_6/appearance_metric_sensitivity.csv').exists():
        return pd.read_csv(OUT/'tables/stimulus_descriptors.csv')
    a=pd.read_csv(REPO/'results/stage_2_6/appearance_metric_sensitivity.csv')
    g=pd.read_csv(REPO/'results/stage_2_6/geometric_metric_reproduced_values.csv')
    return a.drop(columns=['G']).merge(g[['identity','picture_filename','G_new']],on=['identity','picture_filename'],validate='one_to_one')
def behavior():
    if os.environ.get('IPM_ROUND1_REPLAY')=='1' or not (DATA/'ipm_stage2_analysis/02_behavior_primary/behavior_trialwise_long.csv').exists():
        return pd.read_csv(OUT/'tables/behavior_input.csv')
    b=pd.read_csv(DATA/'ipm_stage2_analysis/02_behavior_primary/behavior_trialwise_long.csv')
    b=b.merge(metrics(),left_on=['Identity']+F,right_on=['identity','fslim','eye','mouth','skin'],validate='many_to_one',how='left')
    assert len(b)==23040 and b.A2.notna().all()
    return b.rename(columns={'subj':'Participant','picture_filename':'Stimulus'})
def eeg():
    if os.environ.get('IPM_ROUND1_REPLAY')=='1' or not (DATA/'CHB_multimodal_facial_editing/multimodal_trial_master.csv').exists():
        return pd.read_csv(OUT/'tables/eeg_analysis_input.csv')
    rois=[f'{r}_{t}_Mean' for r in ['Occipitotemporal','Posterior','Centroparietal','Frontal'] for t in TAGS]
    cols=['Subject','Picture','Identity','IsOriginal','IsAttentionCheck','EpochAccepted','ArtifactFlag','TrialOrder','Beauty','Naturalness']+F+rois
    d=pd.read_csv(DATA/'CHB_multimodal_facial_editing/multimodal_trial_master.csv',usecols=cols,low_memory=False)
    d=d[d.Subject.isin(SUB_EEG)&d.IsOriginal.eq(0)&d.IsAttentionCheck.eq(0)&d.EpochAccepted.eq(True)&d.ArtifactFlag.eq(False)].copy()
    im=metrics()
    for c in M:im[c]=(im[c]-im[c].mean())/im[c].std(ddof=0)
    d=d.merge(im[['identity','picture_filename']+M+COMP],left_on=['Identity','Picture'],right_on=['identity','picture_filename'],validate='many_to_one',how='left')
    d=d.dropna(subset=['TrialOrder','Beauty','Naturalness']+F+M+rois)
    assert len(d)==9626 and d.Subject.nunique()==28
    return d
def formula(extra,identity=True):
    return 'Rating ~ RatingDimension_Beauty*('+'+'.join(F+extra)+')'+('+C(Identity)' if identity else '')
def fit_predict(train,test,extra,identity=True):
    y,x=patsy.dmatrices(formula(extra,identity),train,return_type='dataframe')
    xt=patsy.build_design_matrices([x.design_info],test)[0]
    assert np.linalg.matrix_rank(x)==x.shape[1]
    return np.asarray(xt)@np.linalg.lstsq(x,np.asarray(y).ravel(),rcond=None)[0]
def scores(y,p,w=None):
    y=np.asarray(y);p=np.asarray(p)
    avg=lambda v:np.average(v,weights=w)
    return {'r2':1-avg((y-p)**2)/avg((y-avg(y))**2),'rmse':np.sqrt(avg((y-p)**2)),'mae':avg(np.abs(y-p))}
def folds():
    s=np.array(SUB);np.random.default_rng(20260818).shuffle(s)
    return np.array_split(s,5)
def eeg_betas(d,descriptor=True,component=False):
    rows=[];diagnostics=[]
    terms=F+(['G_new']+COMP+['I'] if component else M) if descriptor else F
    for s in SUB_EEG:
        q=d[d.Subject.eq(s)].copy()
        order=q.TrialOrder;z=(order-order.mean())/order.std(ddof=0)
        x=pd.concat([pd.Series(1.,index=q.index,name='Intercept'),q[terms].astype(float),z.rename('TrialOrder_z'),pd.get_dummies(q.Identity,prefix='id',drop_first=True,dtype=float)],axis=1)
        a=x.to_numpy();assert np.linalg.matrix_rank(a)==a.shape[1]
        b=np.linalg.lstsq(a,q[EC].to_numpy(),rcond=None)[0]
        for j,w in enumerate(W):
            for k,t in enumerate(x.columns):
                if t in (['G_new','I']+COMP if component else M) if descriptor else t in F:
                    rows.append({'Subject':s,'window':w,'metric':t,'beta_uV':b[k,j]})
        std=a[:,1:].std(axis=0);az=(a[:,1:]-a[:,1:].mean(axis=0))/std
        hat=np.sum(a*np.linalg.pinv(a).T,axis=1)
        diagnostics.append({'Subject':s,'n':len(q),'p':a.shape[1],'condition_scaled':np.linalg.cond(az),'max_leverage':hat.max()})
    return pd.DataFrame(rows),pd.DataFrame(diagnostics)
def group(b,seed=20260904):
    subjects=[s for s in SUB_EEG if s in set(b.Subject)]
    labels=[];vs=[]
    for key,q in b.groupby(['window','metric'],sort=False):
        labels.append(key);vs.append(q.set_index('Subject').reindex(subjects).beta_uV.to_numpy())
    a=np.column_stack(vs);n=len(a);t=stats.ttest_1samp(a,0,axis=0).statistic
    rng=np.random.default_rng(seed);null=[]
    for start in range(0,10000,500):
        v=a[None]*rng.choice([-1.,1.],size=(500,n,1))
        tv=v.mean(1)/(v.std(1,ddof=1)/np.sqrt(n));null.extend(np.max(np.abs(tv),axis=1))
    null=np.array(null);rows=[]
    for j,(w,m) in enumerate(labels):
        v=a[:,j];se=v.std(ddof=1)/np.sqrt(n);mu=v.mean();crit=stats.t.ppf(.975,n-1)
        rows.append(dict(window=w,metric=m,n_participants=n,mean_beta=mu,se=se,ci95_low=mu-crit*se,ci95_high=mu+crit*se,t=t[j],p_uncorrected=2*stats.t.sf(abs(t[j]),n-1),p_maxT_familywise=(1+np.sum(null>=abs(t[j])))/10001,cohen_dz=mu/v.std(ddof=1),family_tests=a.shape[1],seed=seed))
    return pd.DataFrame(rows)
