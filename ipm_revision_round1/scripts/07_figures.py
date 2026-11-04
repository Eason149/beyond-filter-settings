from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.titleweight':'bold'})
coral='#D66C61';blue='#477E9C';gold='#C59B48';gray='#66717B'
labels=['Operations','+ Luminance','+ Chromatic','+ Texture','+ A2','+ Three components']
def export(fig,name):
    fig.savefig(OUT/'figures'/f'{name}.png',dpi=250,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/'figures'/f'{name}.svg',bbox_inches='tight',facecolor='white');plt.close(fig)
a=pd.read_csv(OUT/'a2_ablation_results.csv');a=a[a.dimension.eq('pooled')]
u=pd.read_csv(OUT/'tables/ablation_paired_uncertainty.csv')
fig,axs=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
for j,m in enumerate(MODELS):
    q=a[a.model.eq(m)&a.fold.ne('mean')];mu=a[a.model.eq(m)&a.fold.eq('mean')].r2.iloc[0]
    axs[0].scatter(q.r2,np.full(len(q),j)+np.linspace(-.10,.10,len(q)),color=gray,s=17,alpha=.55)
    axs[0].scatter(mu,j,s=62,color=coral if m=='M4' else blue,zorder=3)
    axs[0].text(mu+.003,j+.23,f'{mu:.4f}',fontsize=9,color=gray)
axs[0].set(yticks=np.arange(6),yticklabels=labels,xlabel='Held-out-participant R² (pooled ratings)',title='A   Ablation: same five folds');axs[0].invert_yaxis();axs[0].grid(axis='x',alpha=.15)
units=['participant','stimulus_within_identity','participant_x_stimulus'];names=['Participants','Images within identity','Participants × images']
for j,unit in enumerate(units):
    r=u[(u.bootstrap_unit==unit)&(u.comparison=='M4 minus M0')&(u.metric=='r2')].iloc[0]
    axs[1].errorbar(r.estimate,j,xerr=[[r.estimate-r.ci95_low],[r.ci95_high-r.estimate]],fmt='o',color=coral,capsize=4)
axs[1].axvline(0,color=gray,lw=1,ls='--');axs[1].set(yticks=range(3),yticklabels=names,xlabel='A2 increment in R²; pointwise 95% interval',title='B   Conditional paired uncertainty');axs[1].invert_yaxis();axs[1].grid(axis='x',alpha=.15)
fig.suptitle('A2 adds information to operations; its advantage depends on the comparison',fontsize=13)
export(fig,'S1_a2_ablation')

im=metrics();ss=pd.read_csv(OUT/'skin_sensitivity_results.csv');dist=pd.read_csv(OUT/'tables/skin_descriptor_distribution.csv')
fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained');colors=[coral,blue,gold,'#64856B']
for k,id in enumerate(IDS):
    for skin in [0,1]:
        q=im[(im.identity==id)&(im.skin==skin)];axs[0].scatter(np.full(len(q),skin)+(k-1.5)*.08+np.linspace(-.025,.025,len(q)),q.A2,s=24,color=colors[k],alpha=.85,label=id if skin==0 else None)
axs[0].set(xticks=[0,1],xticklabels=['Skin = 0','Skin = 1'],ylabel='A2 (global stimulus SD)',title='A   All 64 images');axs[0].legend(frameon=False,ncol=2,fontsize=9);axs[0].grid(axis='y',alpha=.15)
ylabels=[]
for i,(skin,dimension) in enumerate([(0,'unnaturalness'),(0,'beauty'),(1,'unnaturalness'),(1,'beauty')]):
    q=ss[(ss.model==f'Skin{skin}')&(ss.dimension==dimension)].iloc[0]
    sd=dist[(dist.Skin==skin)&(dist.identity=='all')&(dist.descriptor=='A2')].sd.iloc[0]
    est=q.estimate*sd;lo=q.ci95_low*sd;hi=q.ci95_high*sd
    axs[1].errorbar(est,i,xerr=[[est-lo],[hi-est]],fmt='o',color=coral if dimension=='beauty' else blue,capsize=4)
    ylabels.append(f'Skin {skin} · {dimension}')
axs[1].axvline(0,color=gray,lw=1,ls='--');axs[1].set(yticks=range(4),yticklabels=ylabels,xlabel='Rating difference per within-Skin SD of A2',title='B   Adjusted within-Skin associations');axs[1].invert_yaxis();axs[1].grid(axis='x',alpha=.15)
fig.suptitle('Within-Skin variation carries information with different association directions',fontsize=13)
export(fig,'S2_skin_sensitivity')

h=pd.read_csv(OUT/'identity_holdout_results.csv');h=h[h.dimension.eq('pooled')]
fig,axs=plt.subplots(1,2,figsize=(11.5,5.2),layout='constrained')
for ax,value,title in [(axs[0],'r2','A   Held-out-identity R²'),(axs[1],'delta_r2_vs_M0','B   Difference from operations')]:
    grid=h.pivot(index='model',columns='held_out_identity',values=value).reindex(index=list(MODELS),columns=IDS).to_numpy();v=np.max(abs(grid));img=ax.imshow(grid,cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto')
    for i in range(6):
        for j in range(4):ax.text(j,i,f'{grid[i,j]:+.3f}',ha='center',va='center',fontsize=10,color='white' if abs(grid[i,j])>v*.65 else '#1F2830')
    ax.set(xticks=range(4),xticklabels=IDS,yticks=range(6),yticklabels=labels,title=title,xlabel='Identity excluded from all training outcomes')
    fig.colorbar(img,ax=ax,shrink=.8,pad=.03)
fig.suptitle('Transfer across four source identities: each fold retained',fontsize=13)
export(fig,'S3_identity_holdout')

e=pd.read_csv(OUT/'eeg_identity_robustness.csv');e=e[e.metric.eq('I')]
fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True,layout='constrained')
for ax,w in zip(axs,W):
    selections=[('All identities','baseline','none')]+[('Omit '+id,'leave_identity_out',id) for id in IDS]+[('Components replace A2','components_replace_A2','A2 replaced by three components')]
    labs=[]
    for j,(lab,analysis,omitted) in enumerate(selections):
        r=e[(e.window==w)&(e.analysis==analysis)&(e.omitted==omitted)].iloc[0];color=coral if analysis=='baseline' else blue
        ax.errorbar(r.mean_beta,j,xerr=[[r.mean_beta-r.ci95_low],[r.ci95_high-r.mean_beta]],fmt='o',color=color,capsize=4)
        labs.append(lab);ax.text(.42,j,f'pFWE={r.p_maxT_familywise:.4f}',va='center',fontsize=8,color=gray)
    ax.axvline(0,color=gray,lw=1,ls='--');ax.set(xlim=(-.65,.68),yticks=range(len(labs)),yticklabels=labs,title=w,xlabel='I coefficient (µV per global descriptor SD)');ax.grid(axis='x',alpha=.15)
axs[0].invert_yaxis();fig.suptitle('I–EEG association: identity omission and component adjustment',fontsize=13)
export(fig,'S4_eeg_identity_sensitivity')

c=im[['G_new','A2','I']+COMP+['fslim','eye','mouth','skin']].corr()
labs=['G','A2','I','Luminance','Chromatic','Texture','Slimming','Eye','Mouth','Skin']
fig,ax=plt.subplots(figsize=(9,8),layout='constrained');img=ax.imshow(c,cmap='RdBu_r',vmin=-1,vmax=1)
for i in range(len(c)):
    for j in range(len(c)):
        ax.text(j,i,f'{c.iloc[i,j]:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(c.iloc[i,j])>.65 else '#1F2830')
ax.set(xticks=range(10),xticklabels=labs,yticks=range(10),yticklabels=labs,title='Image-level relationships among descriptors, components and operations (n = 64)');plt.setp(ax.get_xticklabels(),rotation=45,ha='right');fig.colorbar(img,ax=ax,shrink=.8,label='Pearson r')
export(fig,'S5_image_correlations')
(OUT/'figures/README.md').write_text('''# Round 1 supplementary figures

All figures are generated from completed Round 1 tables by `scripts/07_figures.py`. PNG is 250 dpi; SVG preserves editable vector labels. These files do not replace manuscript artwork.

- S1: Left points are the five original participant folds, large dots their mean. Right intervals are paired bootstrap evaluation-sample uncertainty conditional on fixed OOF predictions; participant, identity-stratified image, and product participant/image resampling. They do not include refitting uncertainty or population uncertainty across identities. Intervals are pointwise, not simultaneous.
- S2: Each dot represents one edited image. The adjusted A2 slopes come from crossed participant and stimulus intercept/dimension-slope models within each Skin group. Plot slopes and Wald 95% intervals are converted from global A2 units to the observed within-group SD for interpretability; original units and Holm-adjusted p-values are retained in `skin_sensitivity_results.csv`. Association direction differs across Skin groups; no causal interpretation.
- S3: Every identity-held-out fold; outcome-independent A2/component transformations are fitted using training images only. Both outcomes are pooled. Negative R² means worse than the mean outcome in the test fold; no identity dummies or held-out outcomes enter training. This is identity transfer with participant overlap.
- S4: Point estimates and unadjusted t-based 95% intervals for participant-level I coefficients. Baseline/identity-omission p-values are six-test max-|t| corrected (G,A2,I × two windows); the exploratory component-replacement p-values use a two-test I family. The figure displays different prespecified families explicitly in this caption. Coefficient scales and predefined ROI/windows are unchanged.
- S5: Pearson correlations across 64 stimuli, not trial-level correlation counts. Spearman and partial correlations are in the accompanying tables. This descriptive heatmap has no significance annotations.
''',encoding='utf-8')
print('Five supplementary figures exported as PNG and SVG.')
