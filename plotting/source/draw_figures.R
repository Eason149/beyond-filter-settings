# Compact scientific figures from archived results; no model refitting.
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep('^--file=', args, value=TRUE)
script <- normalizePath(sub('^--file=', '', file_arg[1]), winslash='/')
root <- Sys.getenv('IPM_FIGURE_ROOT', unset=dirname(dirname(script)))
for (subdir in c('panels','quality_checks','review_layouts')) dir.create(file.path(root,subdir),recursive=TRUE,showWarnings=FALSE)
local_lib <- Sys.getenv('FIGURE_R_LIBRARY', unset='')
if (nzchar(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(grid); library(jsonlite)})
skill <- Sys.getenv('NATURE_FIGURE_SCRIPTS', unset=file.path(dirname(dirname(script)),'tools'))
source(file.path(dirname(script),'panel_alignment.R'))
python <- Sys.getenv('FIGURE_QA_PYTHON', unset=Sys.which('python'))
data_dir <- file.path(root,'data'); qa_dir <- file.path(root,'quality_checks')
read_data <- function(name) read.csv(file.path(data_dir,name),check.names=FALSE,stringsAsFactors=FALSE)
CORAL <- '#C56655'; BLUE <- '#42788F'; TEAL <- '#5C9B8C'; LILAC <- '#9382AB'; INK <- '#29323A'
pal <- c(FSlim=CORAL,Eye=BLUE,Mouth=LILAC,Skin=TEAL)
metric_pal <- c(G_new=TEAL,A2=CORAL,I=BLUE)
ops <- c('FSlim','Eye','Mouth','Skin'); op_labels <- c('Slimming','Eye size','Mouth size','Skin')
metrics <- c('G_new','A2','I'); windows <- c('350-600 ms','600-1000 ms')
theme_set(theme_classic(base_size=7,base_family='sans') + theme(
 axis.title=element_text(size=6.7,colour=INK),axis.text=element_text(size=6.3,colour=INK),
 axis.line=element_line(linewidth=.25,colour=INK),axis.ticks=element_line(linewidth=.22),
 axis.ticks.length=unit(1,'mm'),axis.title.x=element_text(margin=margin(t=2)),
 axis.title.y=element_text(margin=margin(r=2)),
 plot.title=element_text(size=7.6,face='bold',colour=INK,margin=margin(b=4)),
 plot.tag=element_text(size=9,face='bold'),plot.tag.position='topleft',
 plot.margin=margin(3,4,3,3,unit='pt'),
 legend.title=element_text(size=6.1),legend.text=element_text(size=5.8),
 legend.key.height=unit(3,'mm'),legend.key.width=unit(2.4,'mm'),
 legend.margin=margin(0,0,0,0),legend.box.spacing=unit(1,'mm'),legend.background=element_blank(),
 panel.grid=element_blank()))
pf <- function(x) ifelse(x<.001,'< .001',sub('^0','',sprintf('%.3f',x)))
plab <- function(x) ifelse(x<.001,'p < .001',paste0('p = ',pf(x)))
zero <- function() geom_hline(yintercept=0,colour='#B9BEC4',linetype='dashed',linewidth=.24)
tag <- function(p,i) p+labs(tag=letters[i])
stim <- read_data('stimulus_plot_data.csv'); stopifnot(nrow(stim)==64)
ob <- read_data('behavior_operation_subject_betas.csv'); oe <- read_data('predefined_window_subject_betas.csv')
ogb <- read_data('behavior_operation_group_stats.csv'); oge <- read_data('predefined_window_group_stats.csv')
gg <- read_data('gai_to_eeg_group_stats.csv'); subj <- read_data('participant_GAI_EEG_coefficients.csv')
base <- read_data('eeg_to_behavior_simple_subject_betas.csv'); adj <- read_data('eeg_to_behavior_adjusted_subject_betas.csv')
gb <- read_data('eeg_to_behavior_simple_group_stats.csv'); ga <- read_data('eeg_to_behavior_adjusted_group_stats.csv')
stopifnot(nrow(subj)==28)

dist_plot <- function(raw,stats,levels,labels,colours,title,ylab,limits){
 raw$group <- factor(raw$group,levels=levels); stats$group <- factor(stats$group,levels=levels)
 stopifnot(all(table(raw$group)==28),all(is.finite(raw$value)))
 means <- tapply(raw$value,raw$group,mean)
 stopifnot(all(abs(means[as.character(stats$group)]-stats$mean_beta)<1e-10))
 stats$fillpoint <- ifelse(stats$p_maxT_familywise<.05,colours[as.character(stats$group)],'white')
 span <- diff(limits); stats$py <- limits[2]-.035*span
 ggplot(raw,aes(group,value,fill=group)) + zero() +
  geom_violin(width=.64,trim=TRUE,bw='nrd0',alpha=.18,colour=NA) +
  geom_point(aes(colour=group),position=position_jitter(width=.12,height=0,seed=37),size=.80,alpha=.66,stroke=0) +
  geom_errorbar(data=stats,aes(x=group,ymin=ci95_low,ymax=ci95_high),inherit.aes=FALSE,width=.12,linewidth=.42,colour=INK) +
  geom_point(data=stats,aes(x=group,y=mean_beta),inherit.aes=FALSE,shape=21,size=2.05,stroke=.48,fill=stats$fillpoint,colour=INK) +
  geom_text(data=stats,aes(x=group,y=py,label=plab(p_maxT_familywise)),inherit.aes=FALSE,size=2.05,colour=INK) +
  scale_fill_manual(values=colours,guide='none')+scale_colour_manual(values=colours,guide='none')+
  scale_x_discrete(labels=labels)+scale_y_continuous(limits=limits,expand=expansion(mult=0))+
  labs(x=NULL,y=ylab,title=title)
}

# Figure 4: complementary measures and the full stimulus design.
R <- cor(stim[,c('G_new','A2','I')]); cr <- expand.grid(x=1:3,y=1:3)
cr$r <- as.vector(R);cr <- cr[cr$x<=cr$y,]
p4a <- ggplot(cr,aes(x,y,fill=r))+geom_tile()+
 geom_text(aes(label=sprintf('%.2f',r)),size=2.5,colour=ifelse(abs(cr$r)>.7,'white',INK))+
 scale_fill_gradient2(low=BLUE,mid='#FCFCFA',high=CORAL,limits=c(-1,1),name='Pearson r',breaks=c(-1,0,1))+
 scale_x_continuous(breaks=1:3,labels=c('G','A2','I'),expand=c(0,0))+
 scale_y_reverse(breaks=1:3,labels=c('G','A2','I'),expand=c(0,0))+
 guides(fill=guide_colorbar(barheight=unit(18,'mm'),barwidth=unit(2.4,'mm'),display='rectangles'))+
 labs(x=NULL,y=NULL,title='Relationships among G, A2 and I')+theme(axis.line=element_blank(),axis.ticks=element_blank())
stim$skin_group <- factor(stim$skin,levels=c(0,1),labels=c('Skin 0','Skin 1'))
p4b <- ggplot(stim,aes(A,A2,colour=skin_group,shape=identity))+
 geom_point(size=1.55,alpha=.88,stroke=.25)+
 scale_colour_manual(values=c(BLUE,CORAL),name=NULL)+scale_shape_manual(values=c(16,15,17,18),guide='none')+
 annotate('text',x=-1.45,y=2.28,label=sprintf('r = %.3f',cor(stim$A,stim$A2)),hjust=0,size=2.2)+
 scale_x_continuous(breaks=c(-1,0,1,2),limits=c(-1.55,2.7))+
 scale_y_continuous(breaks=c(-1,0,1,2),limits=c(-1.2,2.5))+
 labs(x='Broad appearance A (z)',y='Local appearance A2 (z)',title='Broad and local appearance')+
 theme(legend.position='inside',legend.position.inside=c(.8,.9),legend.background=element_blank(),legend.key.height=unit(2.5,'mm'))
design <- do.call(rbind,lapply(seq_along(ops),function(j)data.frame(image=1:64,row=j,value=stim[[c('fslim','eye','mouth','skin')[j]]])))
p4c <- ggplot(design,aes(image,row,fill=factor(value)))+geom_tile()+
 geom_vline(xintercept=c(16.5,32.5,48.5),colour='white',linewidth=.6)+
 scale_fill_manual(values=c('#F1F2F3','#626D77'),guide='none')+
 scale_x_continuous(breaks=c(8.5,24.5,40.5,56.5),labels=c('F1','F2','M1','M2'),expand=c(0,0))+
 scale_y_reverse(breaks=1:4,labels=op_labels,expand=c(0,0))+
 labs(x='Source identity',y=NULL,title='Editing conditions · 64 images')+
 theme(axis.line=element_blank(),axis.ticks=element_blank())
Z <- scale(stim[,metrics],center=TRUE,scale=apply(stim[,metrics],2,function(x)sqrt(mean((x-mean(x))^2))))
zm <- do.call(rbind,lapply(1:3,function(j)data.frame(image=1:64,row=j,value=Z[,j])))
lim <- ceiling(max(abs(zm$value)))
p4d <- ggplot(zm,aes(image,row,fill=value))+geom_tile()+
 geom_vline(xintercept=c(16.5,32.5,48.5),colour='white',linewidth=.6)+
 scale_fill_gradient2(low='#28647C',mid='#FAFAF7',high='#BA5748',limits=c(-lim,lim),name='z score',breaks=c(-lim,0,lim))+
 scale_x_continuous(breaks=c(8.5,24.5,40.5,56.5),labels=c('F1','F2','M1','M2'),expand=c(0,0))+
 scale_y_reverse(breaks=1:3,labels=c('G','A2','I'),expand=c(0,0))+
 guides(fill=guide_colorbar(barheight=unit(12,'mm'),barwidth=unit(2.4,'mm'),display='rectangles'))+
 labs(x='Source identity',y=NULL,title='Realized visual changes')+
 theme(axis.line=element_blank(),axis.ticks=element_blank())

# Figure 5: acquisition, stimulus presentation and temporal analysis.
cap <- read_data('acquisition_64_channels.csv');stopifnot(nrow(cap)==64,sum(cap$selected)==6)
th <- seq(0,2*pi,length.out=361);head <- data.frame(x=1.24*cos(th),y=1.24*sin(th))
nose <- data.frame(x=c(-.12,0,.12),y=c(1.235,1.38,1.235))
ears <- do.call(rbind,lapply(c(-1,1),function(s)data.frame(x=s*(1.24+.065*sin(seq(0,pi,length.out=60))),y=.17*cos(seq(0,pi,length.out=60)),side=s)))
cap$tx <- ifelse(cap$channel%in%c('M1','M2'),sign(cap$x)*1.43,cap$x)
cap$col <- ifelse(cap$selected,CORAL,'#6B7580')
p5a <- ggplot()+geom_path(data=head,aes(x,y),colour='#8B929A',linewidth=.35)+
 geom_path(data=nose,aes(x,y),colour='#8B929A',linewidth=.35)+geom_path(data=ears,aes(x,y,group=side),colour='#8B929A',linewidth=.35)+
 geom_point(data=cap,aes(x,y),colour=cap$col,size=ifelse(cap$selected,2.2,.8))+
 geom_text(data=cap,aes(tx,y+.085,label=channel),size=1.95,colour=cap$col,vjust=0)+
 coord_cartesian(xlim=c(-1.62,1.62),ylim=c(-1.43,1.55),expand=FALSE)+
 labs(title='64-channel montage')+theme_void(base_family='sans')+
 theme(plot.title=element_text(size=7.6,face='bold',margin=margin(b=4)),plot.margin=margin(3,4,3,3),plot.tag=element_text(size=9,face='bold'))
timeline <- function(start,end,labs,cols,title,xlim,breaks){
 dat <- data.frame(start=start,end=end,y=rev(seq_along(start)),label=labs,col=cols)
 ggplot(dat)+geom_rect(aes(xmin=start,xmax=end,ymin=y-.17,ymax=y+.17),fill=cols,colour=NA)+
  geom_text(aes(x=(start+end)/2,y=y+.38,label=paste0(start,'\u2013',end)),size=2.1)+
  scale_y_continuous(breaks=dat$y,labels=labs,limits=c(.5,length(start)+.7),expand=c(0,0))+
  scale_x_continuous(limits=xlim,breaks=breaks,expand=c(0,0))+
  labs(x='Time from face onset (ms)',y=NULL,title=title)+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
}
p5b <- timeline(c(-200,0),c(0,1400),c('Baseline','Face'),c('#A5AFB8','#D6DCE1'),'Stimulus presentation',c(-260,1480),c(-200,0,600,1400))
p5c <- timeline(c(350,600,300,500),c(600,1000,500,1000),c('Study: middle','Study: late','Ref: P300','Ref: LPP'),c(BLUE,CORAL,'#BCD2DD','#E8BEB6'),'Analysis windows',c(220,1080),c(300,600,1000))

# Figure 6: distributions of all 28 archived participant operation effects.
p6 <- lapply(1:4,function(k){
 if(k<=2){o<-c('Beauty','Naturalness')[k];v<-ob[ob$outcome==o,];s<-ogb[ogb$outcome==o,];v$value<-v$rating_beta;ylab<-'Rating difference';limits<-c(-1.9,3.2)}
 else{w<-c('MiddleLate_350_600','Late_600_1000')[k-2];v<-oe[oe$window==w,];s<-oge[oge$window==w,];v$value<-v$beta_uV;ylab<-'Amplitude difference (\u00b5V)';limits<-c(-1.55,1.35)}
 v$group<-v$factor;s$group<-s$factor
 dist_plot(v,s,ops,op_labels,pal,c('Beauty','Perceived unnaturalness','EEG: 350\u2013600 ms','EEG: 600\u20131000 ms')[k],ylab,limits)
})

# Figure 8: group distributions above individual coefficient maps.
p8dist <- lapply(1:2,function(k){
 v<-do.call(rbind,lapply(1:3,function(j)data.frame(group=metrics[j],value=subj[[1+(k-1)*3+j]])))
 s<-gg[gg$window==windows[k],];s$group<-s$metric
 dist_plot(v,s,metrics,c('G','A2','I'),metric_pal,gsub('-','\u2013',windows[k]),'EEG coefficient (\u00b5V / SD)',c(-1.5,1.8))
})
p8heat <- lapply(1:2,function(k){
 h<-do.call(rbind,lapply(1:3,function(j)data.frame(metric=j,subject=1:28,beta=subj[[1+(k-1)*3+j]])))
 ggplot(h,aes(metric,subject,fill=beta))+geom_tile()+
 scale_fill_gradient2(low='#28647C',mid='#FAFAF7',high='#BA5748',limits=c(-1.5,1.5),name='\u00b5V / SD',breaks=c(-1.5,0,1.5))+
 scale_y_reverse(breaks=1:28,labels=subj$participant,expand=c(0,0))+
 scale_x_continuous(breaks=1:3,labels=c('G','A2','I'),expand=c(0,0))+
 guides(fill=guide_colorbar(barheight=unit(35,'mm'),barwidth=unit(2.4,'mm'),display='rectangles'))+
 labs(x=NULL,y=NULL,title=paste('Participant profiles:',gsub('-','\u2013',windows[k])))+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=5.3))
})

# Figure 9: archived joint-model estimates and paired participant folds.
bt<-read_data('behavior_frozen_metric_models.csv');bt<-bt[bt$model=='A2',]
p9coef<-lapply(1:2,function(k){
 terms<-paste0(if(k==1)'' else 'RatingDimension_Beauty:',metrics)
 v<-bt[match(terms,bt$term),];stopifnot(!anyNA(v$estimate));v$group<-factor(metrics,levels=rev(metrics));v$lo<-v$estimate-1.96*v$se;v$hi<-v$estimate+1.96*v$se
 ggplot(v,aes(estimate,group,colour=group))+geom_vline(xintercept=0,colour='#B9BEC4',linetype='dashed',linewidth=.24)+
 geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.16,linewidth=.5)+geom_point(size=2.1)+
 scale_colour_manual(values=metric_pal,guide='none')+scale_y_discrete(labels=c(I='I',A2='A2',G_new='G'))+
 scale_x_continuous(limits=c(-.50,.33),breaks=c(-.4,-.2,0,.2))+
 labs(x=if(k==1)'Unnaturalness slope (95% CI)' else 'Beauty minus unnaturalness slope (95% CI)',y=NULL,title=c('Descriptor associations','Differences between rating dimensions')[k])
})
cv<-read_data('a2_behavior_participant_cv.csv');mods<-c('A2_only','factor','factor_plus_A2');cv$x<-match(cv$model,mods)
means<-aggregate(r2~x,cv,mean)
p9c<-ggplot(cv,aes(x,r2,group=fold))+geom_line(colour='#BAC1C7',linewidth=.3)+geom_point(shape=21,size=1.35,fill='white',colour='#9CA7B0',stroke=.3)+
 geom_line(data=means,aes(x,r2,group=1),inherit.aes=FALSE,colour=CORAL,linewidth=.7)+geom_point(data=means,aes(x,r2),inherit.aes=FALSE,shape=15,colour=CORAL,size=2)+
 geom_text(data=means,aes(x,y=-.035,label=sprintf('%.4f',r2)),inherit.aes=FALSE,size=2.1,colour=CORAL,fontface='bold')+
 scale_x_continuous(breaks=1:3,labels=c('A2 only','Operations','Operations\n+ A2'),limits=c(.7,3.3))+
 scale_y_continuous(limits=c(-.06,.29),breaks=c(0,.1,.2),expand=c(0,0))+
 labs(x=NULL,y='Held-out R\u00b2',title='Participant-held-out prediction')
delta<-merge(cv[cv$model=='factor_plus_A2',c('fold','r2')],cv[cv$model=='factor',c('fold','r2')],by='fold',suffixes=c('_new','_base'));delta$gain<-delta$r2_new-delta$r2_base
delta$row<-6-delta$fold
gain_summary<-data.frame(gain=mean(delta$gain),row=0)
p9d<-ggplot(delta,aes(gain,row))+geom_segment(aes(x=0,xend=gain,yend=row),colour='#D0D6DB',linewidth=.5)+
 geom_point(size=1.7,colour=BLUE)+geom_text(aes(x=gain+.0012,label=sprintf('%.4f',gain)),hjust=0,size=2.1)+
 geom_hline(yintercept=.5,colour='#E3E6E8',linewidth=.25)+
 geom_point(data=gain_summary,aes(gain,row),inherit.aes=FALSE,shape=18,size=2.9,colour=CORAL)+
 geom_text(data=gain_summary,aes(x=gain+.0012,y=row,label=sprintf('+%.4f',gain)),inherit.aes=FALSE,hjust=0,size=2.2,colour=CORAL,fontface='bold')+
 scale_x_continuous(limits=c(0,.034),breaks=c(0,.01,.02,.03),expand=expansion(mult=c(.02,.01)))+
 scale_y_continuous(breaks=0:5,labels=c('Mean',paste('Fold',5:1)),expand=expansion(add=.45))+
 labs(x='\u0394R\u00b2: adding A2 to operations',y=NULL,title='A2 improves prediction in every fold')

# Figure 10: four panels, preserving every base/adjusted participant pair.
rating_col<-'rating_beta_per_within_subject_eeg_sd'
p10<-lapply(1:4,function(k){
 o<-if(k<=2)'Beauty' else 'Naturalness';j<-if(k%%2==1)1 else 2;w<-windows[j];co<-c(BLUE,CORAL)[j]
 b<-base[base$outcome==o & base$window==w,];a<-adj[adj$outcome==o & adj$window==w,]
 stopifnot(nrow(b)==28,nrow(a)==28,setequal(b$Subject,a$Subject))
 b$model_x<-1;a$model_x<-2;v<-rbind(b,a);v$value<-v[[rating_col]]
 s1<-gb[gb$outcome==o & gb$window==w,];s2<-ga[ga$outcome==o & ga$window==w,];s<-rbind(s1,s2);s$model_x<-1:2
 stopifnot(abs(mean(b[[rating_col]])-s1$mean_beta)<1e-10,abs(mean(a[[rating_col]])-s2$mean_beta)<1e-10)
 s$fc<-ifelse(s$p_maxT_familywise<.05,co,'white')
 ggplot(v,aes(model_x,value,group=Subject))+zero()+geom_line(colour='#C7CDD2',alpha=.60,linewidth=.22)+
 geom_point(colour='#8997A1',size=.85,alpha=.50)+
 geom_line(data=s,aes(x=model_x,y=mean_beta,group=1),inherit.aes=FALSE,colour=co,linewidth=.65)+
 geom_errorbar(data=s,aes(x=model_x,ymin=ci95_low,ymax=ci95_high),inherit.aes=FALSE,width=.10,linewidth=.5,colour=INK)+
 geom_point(data=s,aes(x=model_x,y=mean_beta),inherit.aes=FALSE,shape=22,size=2.4,stroke=.5,fill=s$fc,colour=INK)+
 geom_text(data=s,aes(x=model_x,y=.265,label=plab(p_maxT_familywise)),inherit.aes=FALSE,size=2.1)+
 scale_x_continuous(breaks=1:2,labels=c('Base','Adjusted'),limits=c(.65,2.35))+
 scale_y_continuous(limits=c(-.24,.29),breaks=c(-.2,0,.2),expand=c(0,0))+
 labs(x=NULL,y='Rating difference / EEG SD',title=paste(if(o=='Beauty')'Beauty' else 'Unnaturalness','\u00b7',gsub('-','\u2013',w)))
})

# Figure 7: unchanged archived channel-time coefficients, shared diverging scale.
channel_time<-read_data('operation_channel_time_means.csv')
channel_order<-unique(channel_time$channel)
stopifnot(nrow(channel_time)==10200,length(channel_order)==50,length(unique(channel_time$time_ms))==51)
stopifnot(!anyDuplicated(channel_time[,c('operation','channel','time_ms')]),all(is.finite(channel_time$mean_beta_uV)))
channel_time$row<-match(channel_time$channel,channel_order)
channel_rows<-c(1,5,9,14,18,23,27,32,36,41,45,50)
map_limit<-ceiling(max(abs(channel_time$mean_beta_uV))*10)/10
stopifnot(abs(map_limit-.7)<1e-12)
# Archived within-operation cluster masks, verified against all original cluster records.
cluster_masks<-read_data('operation_channel_time_cluster_masks.csv')
cluster_masks$row<-match(cluster_masks$channel,channel_order)
stopifnot(nrow(cluster_masks)==10200,sum(cluster_masks$cluster_mask)==252)
p7<-lapply(seq_along(ops),function(k){
 v<-channel_time[channel_time$operation==ops[k],]
 mask<-cluster_masks[cluster_masks$operation==ops[k],]
 contour_layer<-if(any(mask$cluster_mask==1)) geom_contour(data=mask,
  aes(x=time_ms,y=row,z=cluster_mask),breaks=.5,inherit.aes=FALSE,
  colour='black',linewidth=.35,lineend='round',linejoin='round') else NULL
 ggplot(v,aes(time_ms,row,fill=mean_beta_uV))+
 geom_rect(aes(xmin=pmax(0,time_ms-10),xmax=pmin(1000,time_ms+10),ymin=row-.5,ymax=row+.5),colour=NA)+
 geom_vline(xintercept=seq(100,900,100),colour='#39434B',alpha=.14,linewidth=.2)+
 contour_layer+
 scale_fill_gradient2(low='#28647C',mid='#FAFAF7',high='#BA5748',midpoint=0,limits=c(-map_limit,map_limit),breaks=c(-.6,-.3,0,.3,.6),name='Mean coefficient\n(\u00b5V)')+
 scale_x_continuous(breaks=seq(0,1000,200),expand=c(0,0))+
 scale_y_continuous(breaks=channel_rows,labels=gsub('Z','z',channel_order[channel_rows]),expand=c(0,0))+
 coord_cartesian(xlim=c(0,1000),ylim=c(.5,50.5),expand=FALSE)+
 guides(fill=guide_colorbar(barheight=unit(35,'mm'),barwidth=unit(2.4,'mm'),display='rectangles'))+
 labs(x='Time from face onset (ms)',y='Electrode',title=c('Facial slimming','Eye size','Mouth size','Skin editing')[k])+
 theme(axis.line=element_blank(),panel.border=element_rect(fill=NA,colour='#7A848C',linewidth=.25),axis.text.y=element_text(size=5.8))
})

# Evidence roles: stimulus design and measurement; acquisition and timing;
# operations and responses; group EEG estimates and participant variation;
# descriptor associations and held-out gain; paired EEG-rating associations.
panels<-list(`4`=list(p4a,p4b,p4c,p4d),`5`=list(p5a,p5b,p5c),`6`=p6,`7`=p7,`8`=c(p8dist,p8heat),`9`=c(p9coef,list(p9c,p9d)),`10`=p10)
headlines<-c(
 `4`='From editing combinations to realized visual changes',
 `5`='Where and when viewing-period EEG was measured',
 `6`='Editing operations: behavioral and EEG responses',
 `7`='Spatiotemporal editing-operation coefficients',
 `8`='Identity-related displacement is associated with EEG in both windows',
 `9`='A2 adds predictive information beyond editing operations',
 `10`='Late EEG remains associated with beauty after adjustment')
captions<-list(
 `4`='64 images; 16 conditions per identity. Design: 0 light / 1 dark. GAI heatmap: descriptor z scores.',
 `5`='64 channels; coral marks the six-channel ROI. Locations are schematic. Reference windows: Arai & Nittono (2022).',
 `6`='n = 28. All individual coefficients; violin densities; means and 95% CIs. Numbers: familywise p; filled means: p < .05.',
 `7`='Descriptive coefficients: 50 channels, 51 time samples; common scale \u22120.7 to +0.7 \u00b5V.\nBlack contours: exploratory clusters, within-operation corrected p < .05 (Eye and Skin).',
 `8`='n = 28. G: geometry; A2: local appearance; I: identity. Means and 95% CIs; familywise p; shared heatmap color scale.',
 `9`='a\u2013b: original joint model, n = 30; 95% Wald CIs. c\u2013d: five paired folds; coral squares: means; mean \u0394R\u00b2 = 0.0211.',
 `10`='n = 28. Gray: paired participants. Color: means and 95% CIs. Familywise p; filled squares: p < .05.')
heights<-c(`4`=113,`5`=102,`6`=108,`7`=178,`8`=136,`9`=105,`10`=108)
make_layout<-function(num){
 ps<-panels[[num]];ps<-lapply(seq_along(ps),function(i)tag(ps[[i]],i))
 if(num=='7'){
  ps<-lapply(seq_along(ps),function(i){
   p<-ps[[i]]+guides(fill=guide_colorbar(barheight=unit(65,'mm'),barwidth=unit(3,'mm'),display='rectangles'))
   if(i<4) p<-p+labs(x=NULL)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
   p
  })
  p<-wrap_plots(ps,ncol=1,guides='collect')
 }else if(num=='4'){
  # The lower matrices share the full image axis: each column is the same image.
  p<-wrap_plots(ps,design='AB\nCC\nDD',heights=c(1.6,.65,.65))
 }else if(num=='5'){
  # The cap spans the two temporal panels in the right column.
  p<-wrap_plots(ps,design='AB\nAC',widths=c(1.02,1),heights=c(.9,1.1))
 }else{
  p<-wrap_plots(ps,ncol=2,heights=if(num=='8')c(.82,1.45) else if(num=='4')c(1,.68) else c(1,1))
 }
 p+plot_annotation(title=headlines[[num]],caption=captions[[num]],theme=theme(
  plot.title=element_text(size=9.6,face='bold',colour=INK,margin=margin(b=7)),
  plot.caption=element_text(size=5.7,hjust=0,margin=margin(t=3)),plot.margin=margin(3,3,3,3)))
}

save_plot<-function(p,name,width_mm,height_mm,multi=FALSE){
 dest<-file.path(root,if(multi)'review_layouts' else 'panels',name)
 w<-width_mm/25.4;h<-height_mm/25.4
 if(multi){
  require_patchwork_panel_alignment(p,
   manifest_path=file.path(qa_dir,paste0(name,'.alignment-layout.json')),
   report_path=file.path(qa_dir,paste0(name,'.alignment.json')),
   width_in=w,height_in=h,audit_script=file.path(skill,'audit_panel_alignment.py'),python=python,
   tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=TRUE)
 }
 grDevices::cairo_pdf(paste0(dest,'.pdf'),width=w,height=h,family='Arial');print(p);dev.off()
 svglite::svglite(paste0(dest,'.svg'),width=w,height=h);print(p);dev.off()
 ragg::agg_png(paste0(dest,'.png'),width=w,height=h,units='in',res=600);print(p);dev.off()
 ragg::agg_tiff(paste0(dest,'.tif'),width=w,height=h,units='in',res=600,compression='lzw');print(p);dev.off()
 cat('Saved',name,'\n')
}
figures<-list()
requested_figures<-Sys.getenv('FIGURE_ONLY',unset='')
for(num in names(panels)){
 fig<-make_layout(num);figures[[num]]<-fig
 if(nzchar(requested_figures) && !num %in% strsplit(requested_figures,',',fixed=TRUE)[[1]]) next
 save_plot(fig,sprintf('Figure_%02d_review_layout',as.integer(num)),183,heights[[num]],TRUE)
 for(i in seq_along(panels[[num]])){
  h<-if(num=='8' && i>=3)95 else if(num=='5' && i==1)100 else 65
  save_plot(tag(panels[[num]][[i]],i),sprintf('Fig%02d_%s',as.integer(num),toupper(letters[i])),if(num=='7')183 else 92,h)
 }
}
# Preserve each Cairo PDF's exact page dimensions when assembling the atlas.
qpdf::pdf_combine(vapply(names(figures),function(num)file.path(root,'review_layouts',sprintf('Figure_%02d_review_layout.pdf',as.integer(num))),character(1)),output=file.path(root,'Figures_4_to_10_Review.pdf'))
# R-generated contact sheet for review.
ragg::agg_png(file.path(root,'All_figures_preview.png'),width=1800,height=3000,res=300)
grid.newpage();pushViewport(viewport(layout=grid.layout(4,2)))
for(k in seq_along(figures)){
 num<-names(figures)[k];im<-png::readPNG(file.path(root,'review_layouts',sprintf('Figure_%02d_review_layout.png',as.integer(num))))
 aspect<-nrow(im)/ncol(im);w<-min(.98,.95*(3000/4)/(1800/2)/aspect);h<-w*(1800/2)/(3000/4)*aspect
 pushViewport(viewport(layout.pos.row=ceiling(k/2),layout.pos.col=(k-1)%%2+1));grid.raster(im,width=w,height=h,interpolate=TRUE);popViewport()
}
popViewport();dev.off()
cat('COMPLETE: 27 independent panels; 7 compact R figure layouts.\n')
