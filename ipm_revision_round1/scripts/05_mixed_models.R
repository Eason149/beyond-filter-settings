# Round 1 prespecified crossed participant/stimulus dimension-slope models.
# Rscript scripts/05_mixed_models.R <round1 directory> <lme4 library directory>
args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[1],winslash='/',mustWork=TRUE)
if(length(args)>1) .libPaths(c(args[2],.libPaths()))
suppressPackageStartupMessages(library(lme4))
d <- read.csv(file.path(out,'tables/behavior_input.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
d$Participant <- factor(d$Participant);d$Stimulus <- factor(d$Stimulus);d$Identity <- factor(d$Identity)
components <- c('a2_luminance','a2_chroma_js','a2_texture')
u <- d[!duplicated(d$Stimulus),]
for(v in components)d[[v]] <- (d[[v]]-mean(u[[v]]))/sqrt(mean((u[[v]]-mean(u[[v]]))^2))
factors <- c('FSlim','Eye','Mouth','Skin')
extras <- list(M0=character(),M1=components[1],M2=components[2],M3=components[3],M4='A2',M5=components)
random <- '+(1+RatingDimension_Beauty|Participant)+(1+RatingDimension_Beauty|Stimulus)'
makeform <- function(terms)as.formula(paste('Rating ~ RatingDimension_Beauty*(',paste(terms,collapse='+'),')+Identity',random))
infos <- list();coefs <- list();slopes <- list();fits <- list()
runfit <- function(name,dat,terms){
  cat('Fitting',name,'n=',nrow(dat),'\n');flush.console()
  warns <- character()
  fit <- withCallingHandlers(lmer(makeform(terms),data=droplevels(dat),REML=FALSE,
    control=lmerControl(optimizer='bobyqa',calc.derivs=TRUE,optCtrl=list(maxfun=200000,rhobeg=.1,rhoend=1e-8))),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  saveRDS(fit,file.path(out,'models',paste0(name,'.rds')))
  b <- fixef(fit);v <- as.matrix(vcov(fit));se <- sqrt(diag(v));z <- b/se
  coefs[[name]] <<- data.frame(model=name,term=names(b),estimate=unname(b),se=unname(se),ci95_low=unname(b-1.959964*se),ci95_high=unname(b+1.959964*se),p_wald=2*pnorm(-abs(unname(z))))
  msg <- fit@optinfo$conv$lme4$messages
  grad <- tryCatch(max(abs(solve(chol(fit@optinfo$derivs$Hessian),fit@optinfo$derivs$gradient))),error=function(e)NA_real_)
  infos[[name]] <<- data.frame(model=name,n=nobs(fit),n_participants=nlevels(droplevels(dat$Participant)),n_stimuli=nlevels(droplevels(dat$Stimulus)),logLik=as.numeric(logLik(fit)),AIC=AIC(fit),BIC=BIC(fit),singular=isSingular(fit),optimizer_code=fit@optinfo$conv$opt,max_scaled_gradient=grad,convergence=paste(msg,collapse='; '),warnings=paste(warns,collapse='; '))
  if('A2'%in%names(b)){
    for(dim in 0:1){
      contrast <- setNames(rep(0,length(b)),names(b));contrast['A2']<-1
      if(dim==1)contrast['RatingDimension_Beauty:A2']<-1
      est <- sum(contrast*b);s <- sqrt(drop(t(contrast)%*%v%*%contrast))
      slopes[[paste0(name,'_',dim)]] <<- data.frame(model=name,dimension=if(dim==1)'beauty' else 'unnaturalness',estimate=est,se=s,ci95_low=est-1.959964*s,ci95_high=est+1.959964*s,p_wald=2*pnorm(-abs(est/s)))
    }
  }
  write.csv(do.call(rbind,infos),file.path(out,'tables/mixed_fit_diagnostics.csv'),row.names=FALSE)
  write.csv(do.call(rbind,coefs),file.path(out,'tables/mixed_coefficients.csv'),row.names=FALSE)
  if(length(slopes))write.csv(do.call(rbind,slopes),file.path(out,'tables/mixed_A2_slopes.csv'),row.names=FALSE)
  fit
}
for(m in names(extras))fits[[m]] <- runfit(m,d,c(factors,extras[[m]]))
cmp <- function(lhs,rhs){
  a<-fits[[lhs]];b<-fits[[rhs]];lr<-2*(as.numeric(logLik(a))-as.numeric(logLik(b)));df<-attr(logLik(a),'df')-attr(logLik(b),'df')
  data.frame(comparison=paste(lhs,'vs',rhs),LR_chi2=lr,df=df,p=pchisq(lr,df,lower.tail=FALSE),delta_AIC=AIC(a)-AIC(b),delta_BIC=BIC(a)-BIC(b))
}
lr <- do.call(rbind,c(lapply(names(extras)[-1],function(x)cmp(x,'M0')),list(cmp('M5','M4'))));lr$p_holm_six_tests<-p.adjust(lr$p,'holm')
write.csv(lr,file.path(out,'tables/ablation_crossed_comparisons.csv'),row.names=FALSE)
runfit('GAI_full',d,c('G_new','A2','I'))
for(m in c('M0','M4'))fits[[paste0(m,'_valid')]] <- runfit(paste0(m,'_valid'),d[d$Rating>=1 & d$Rating<=7,],c(factors,extras[[m]]))
runfit('GAI_valid',d[d$Rating>=1 & d$Rating<=7,],c('G_new','A2','I'))
write.csv(cmp('M4_valid','M0_valid'),file.path(out,'tables/zero_exclusion_crossed_comparison.csv'),row.names=FALSE)
for(skin in 0:1){
  q<-d[d$Skin==skin,]
  runfit(paste0('Skin',skin),q,c(factors[factors!='Skin'],'A2'))
  for(id in levels(d$Identity))runfit(paste0('Skin',skin,'_omit_',id),q[q$Identity!=id,],c(factors[factors!='Skin'],'A2'))
}
ss<-do.call(rbind,slopes);ss<-ss[grepl('^Skin',ss$model),];ss$p_holm_primary_four<-NA_real_
primary<-ss$model%in%c('Skin0','Skin1');ss$p_holm_primary_four[primary]<-p.adjust(ss$p_wald[primary],'holm')
write.csv(ss,file.path(out,'skin_sensitivity_results.csv'),row.names=FALSE)
extreme<-u$Stimulus[apply(abs(u[,c('G_new','A2','I')])>3,1,any)]
for(i in seq_along(extreme))runfit(paste0('M4_omit_extreme_',i),d[d$Stimulus!=extreme[i],],c(factors,'A2'))
write.csv(data.frame(model=if(length(extreme))paste0('M4_omit_extreme_',seq_along(extreme)) else character(),omitted=extreme),file.path(out,'tables/extreme_model_key.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines('All prespecified mixed-model fits finished.',file.path(out,'logs/mixed_models_complete.txt'))
