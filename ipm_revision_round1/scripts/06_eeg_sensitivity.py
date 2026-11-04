from common import *
import statsmodels.formula.api as smf
d=eeg();im=metrics();allstats=[];allbet=[]
b,diag=eeg_betas(d);g=group(b);g['analysis']='baseline';g['omitted']='none';allstats.append(g)
save(diag,'tables/eeg_design_diagnostics.csv')
for omitted in SUB_EEG:
    g=group(b[~b.Subject.eq(omitted)]);g['analysis']='leave_participant_out';g['omitted']=omitted;allstats.append(g)
print('EEG participant omission completed',flush=True)
for omitted in IDS:
    bet,dg=eeg_betas(d[~d.Identity.eq(omitted)]);g=group(bet);g['analysis']='leave_identity_out';g['omitted']=omitted;allstats.append(g)
    bet['omitted']=omitted;allbet.append(bet)
save(pd.concat(allbet),'tables/eeg_identity_omission_subject_betas.csv')
print('EEG identity omission completed',flush=True)
influence=[]
for j,omitted in enumerate(sorted(d.Picture.unique())):
    bet,_=eeg_betas(d[~d.Picture.eq(omitted)])
    for w in W:
        v=bet.loc[bet.metric.eq('I')&bet.window.eq(w),'beta_uV'];se=v.std(ddof=1)/np.sqrt(len(v));mu=v.mean()
        influence.append(dict(omitted_stimulus=omitted,window=w,n_participants=len(v),mean_beta=mu,se=se,ci95_low=mu-stats.t.ppf(.975,len(v)-1)*se,ci95_high=mu+stats.t.ppf(.975,len(v)-1)*se))
save(pd.DataFrame(influence),'tables/eeg_stimulus_influence.csv')
print('EEG 64-image omission coefficients completed',flush=True)
# Replace A2 by its components, scaled over the same fixed stimulus population.
dc=d.copy()
for c in COMP:dc[c]=(dc[c]-im[c].mean())/im[c].std(ddof=0)
bc,diag=eeg_betas(dc,component=True);save(bc,'tables/eeg_component_adjusted_subject_betas.csv');save(diag,'tables/eeg_component_adjusted_design_diagnostics.csv')
g=group(bc[bc.metric.eq('I')]);g['analysis']='components_replace_A2';g['omitted']='A2 replaced by three components';allstats.append(g)
clean=d[d.Beauty.between(1,7)&d.Naturalness.between(1,7)]
bc,_=eeg_betas(clean);g=group(bc);g['analysis']='valid_ratings_only';g['omitted']=str(len(d)-len(clean))+' zero-rating trials';allstats.append(g)
out=pd.concat(allstats,ignore_index=True);save(out,'eeg_identity_robustness.csv')
im2=im.rename(columns=dict(zip(['fslim','eye','mouth','skin'],F)))
cols=M+COMP+F;corr=[]
for i,a in enumerate(cols):
    for bb in cols[i+1:]:
        corr.append(dict(variable_1=a,variable_2=bb,n_images=len(im),pearson=stats.pearsonr(im2[a],im2[bb]).statistic,spearman=stats.spearmanr(im2[a],im2[bb]).statistic))
save(pd.DataFrame(corr),'tables/stimulus_correlations.csv')
dis=[]
for skin in [0,1]:
    for id in ['all']+IDS:
        q=im2[im2.Skin.eq(skin)];q=q if id=='all' else q[q.identity.eq(id)]
        for c in M+COMP:
            dis.append(dict(Skin=skin,identity=id,descriptor=c,n_images=len(q),mean=q[c].mean(),sd=q[c].std(ddof=1),min=q[c].min(),q25=q[c].quantile(.25),median=q[c].median(),q75=q[c].quantile(.75),max=q[c].max()))
save(pd.DataFrame(dis),'tables/skin_descriptor_distribution.csv')
fit=smf.ols('A2 ~ FSlim+Eye+Mouth+Skin+C(identity)',im2).fit(cov_type='HC3')
save(pd.DataFrame(dict(term=fit.params.index,estimate=fit.params.values,se_HC3=fit.bse.values,p_exploratory=fit.pvalues.values)),'tables/A2_operation_associations.csv')
# Partial image-level correlation after fixed operation and source-identity controls.
resid={c:smf.ols(c+' ~ FSlim+Eye+Mouth+Skin+C(identity)',im2).fit().resid for c in M+COMP}
save(pd.DataFrame([dict(variable=c,partial_r_with_I=stats.pearsonr(resid['I'],resid[c]).statistic,n_images=64) for c in ['G_new','A2']+COMP]),'tables/I_partial_correlations.csv')
(OUT/'logs/eeg_sensitivity_complete.txt').write_text(f'Completed. Baseline n={len(d)}; valid-rating sensitivity n={len(clean)}. Alternative embedding weights unavailable; SFace alone found in bounded search.',encoding='utf-8')
print(out[out.metric.eq('I')].groupby('analysis').agg(beta_min=('mean_beta','min'),beta_max=('mean_beta','max'),p_min=('p_maxT_familywise','min'),p_max=('p_maxT_familywise','max')).to_string(),flush=True)
