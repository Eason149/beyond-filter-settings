# Targeted crossed participant/stimulus robustness analysis for the CHB revision.
# No primary-model replacement, EEG analysis, or cross-validation refit.
args <- commandArgs(trailingOnly=TRUE)
input <- if (length(args)>=1) args[1] else 'results/chb_stimulus_robustness/behavior_stimulus_model_input.csv'
out <- if (length(args)>=2) args[2] else 'results/chb_stimulus_robustness'
if (length(args)>=3) .libPaths(c(args[3],.libPaths()))
suppressPackageStartupMessages(library(lme4))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(input,check.names=FALSE)
d$Participant <- factor(d$Participant)
d$Stimulus <- factor(d$Stimulus)
d$SourceIdentity <- factor(d$SourceIdentity,levels=c('F_1','F_2','M_1','M_2'))
stopifnot(nrow(d)==23040,nlevels(d$Participant)==30,nlevels(d$Stimulus)==64,
          all(complete.cases(d)),all(table(d$Participant,d$Stimulus,d$RatingDimension_Beauty)==6))
forms <- list(
  GAI=Rating ~ RatingDimension_Beauty*(G_new+A2+I)+SourceIdentity+
    (1+RatingDimension_Beauty|Participant)+(1|Stimulus),
  Factor=Rating ~ RatingDimension_Beauty*(FSlim+Eye+Mouth+Skin)+SourceIdentity+
    (1+RatingDimension_Beauty|Participant)+(1|Stimulus),
  Factor_A2=Rating ~ RatingDimension_Beauty*(FSlim+Eye+Mouth+Skin+A2)+SourceIdentity+
    (1+RatingDimension_Beauty|Participant)+(1|Stimulus))
fits <- list(); coefficients <- list(); infos <- list(); variances <- list()
for (name in names(forms)) {
  cat('Fitting',name,'with ML\n');flush.console()
  warns <- character()
  fit <- withCallingHandlers(lmer(forms[[name]],data=d,REML=FALSE,
    control=lmerControl(optimizer='bobyqa',calc.derivs=TRUE,optCtrl=list(maxfun=200000,rhobeg=.1,rhoend=1e-8))),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
  fits[[name]] <- fit;saveRDS(fit,file.path(out,paste0(name,'_ML.rds')))
  b <- fixef(fit);se <- sqrt(diag(vcov(fit)));z <- b/se
  coefficients[[name]] <- data.frame(model=name,term=names(b),estimate=unname(b),
    se=unname(se),ci95_low=unname(b-1.959963984540054*se),ci95_high=unname(b+1.959963984540054*se),
    z=unname(z),p_wald=2*pnorm(-abs(unname(z))),log10_p_wald=(log(2)+pnorm(-abs(unname(z)),log.p=TRUE))/log(10))
  convergence <- fit@optinfo$conv$lme4$messages
  gradient <- fit@optinfo$derivs$gradient
  hessian <- fit@optinfo$derivs$Hessian
  scaled <- tryCatch(max(abs(solve(chol(hessian),gradient))),error=function(e) NA_real_)
  infos[[name]] <- data.frame(model=name,n=nobs(fit),fixed_parameters=length(b),total_parameters=attr(logLik(fit),'df'),
    logLik=as.numeric(logLik(fit)),AIC=AIC(fit),BIC=BIC(fit),sigma=sigma(fit),
    singular=isSingular(fit,tol=1e-4),optimizer_code=fit@optinfo$conv$opt,
    max_abs_gradient=if(is.numeric(gradient)) max(abs(gradient)) else NA_real_,max_abs_scaled_gradient=scaled,
    convergence_message=if(is.null(convergence)) '' else paste(convergence,collapse='; '),
    warnings=paste(warns,collapse='; '))
  v <- as.data.frame(VarCorr(fit));v$model <- name;variances[[name]] <- v
  capture.output(summary(fit),file=file.path(out,paste0(name,'_summary.txt')))
  print(infos[[name]]);flush.console()
}
coef <- do.call(rbind,coefficients);info <- do.call(rbind,infos)
lr <- 2*(as.numeric(logLik(fits$Factor_A2))-as.numeric(logLik(fits$Factor)))
df <- attr(logLik(fits$Factor_A2),'df')-attr(logLik(fits$Factor),'df')
stopifnot(df==2)
comparison <- data.frame(comparison='Factor+A2 vs Factor; crossed stimulus intercept',
  LR_chi2=lr,df=df,p=pchisq(lr,df,lower.tail=FALSE),
  log10_p=pchisq(lr,df,lower.tail=FALSE,log.p=TRUE)/log(10),
  delta_AIC=AIC(fits$Factor_A2)-AIC(fits$Factor),delta_BIC=BIC(fits$Factor_A2)-BIC(fits$Factor))
write.csv(coef,file.path(out,'stimulus_aware_coefficients.csv'),row.names=FALSE)
write.csv(info,file.path(out,'stimulus_aware_fit_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,variances),file.path(out,'stimulus_aware_variance_components.csv'),row.names=FALSE)
write.csv(comparison,file.path(out,'stimulus_aware_A2_increment.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(vapply(forms,function(f)paste(deparse(f),collapse=' '),''),file.path(out,'model_formulas.txt'))
cat('\nNested model comparison:\n');print(comparison)
cat('\nGAI descriptor terms:\n');print(coef[coef$model=='GAI' & grepl('G_new|A2|(^I$)|:I$',coef$term),])
cat('\nIncremental A2 terms:\n');print(coef[coef$model=='Factor_A2' & grepl('A2',coef$term),])
