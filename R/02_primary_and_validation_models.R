project_root <- normalizePath(".", winslash="/", mustWork=TRUE)
out <- file.path(project_root, "results")
base <- file.path(project_root, "data_raw")
dir.create(out, showWarnings=FALSE)
library(survival)
library(survey)
library(splines)
library(ggplot2)
library(foreign)
options(survey.lonely.psu='adjust')
find_xpt <- function(name){z<-list.files(base,pattern=paste0('^',name,'$'),recursive=TRUE,full.names=TRUE,ignore.case=TRUE);stopifnot(length(z)==1);z}
read_keep <- function(name,vars){x<-read.xport(find_xpt(name));x[,intersect(vars,names(x)),drop=FALSE]}
d <- read.csv(file.path(project_root, 'data_processed', 'nhanes_vcf_grip_analysis_dataset.csv'))
demo <- read.xport(find_xpt('DEMO_H.xpt'))[,c('SEQN','DMDEDUC2','INDFMPIR')]
bmx <- read_keep('BMX_H.xpt',c('SEQN','BMXBMI')); smq<-read_keep('SMQ_H.xpt',c('SEQN','SMQ020','SMQ040'))
alq<-read_keep('ALQ_H.xpt',c('SEQN','ALQ101','ALQ110')); paq<-read_keep('PAQ_H.xpt',c('SEQN','PAQ650','PAQ665'))
diq<-read_keep('DIQ_H.xpt',c('SEQN','DIQ010')); bpq<-read_keep('BPQ_H.xpt',c('SEQN','BPQ020'))
mcq<-read_keep('MCQ_H.xpt',c('SEQN','MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F','MCQ220'))
for(x in list(demo,bmx,smq,alq,paq,diq,bpq,mcq)) d<-merge(d,x,by='SEQN',all.x=TRUE)
d$time<-d$PERMTH_EXM/12; d$event<-as.integer(d$MORTSTAT); d$age10<-d$RIDAGEYR/10
d$female<-as.integer(d$RIAGENDR==2); d$race<-factor(d$RIDRETH1); d$education<-factor(ifelse(d$DMDEDUC2 %in% 1:5,d$DMDEDUC2,NA)); d$bmi<-ifelse(d$BMXBMI>10 & d$BMXBMI<100,d$BMXBMI,NA)
d$smoking<-NA_character_; d$smoking[d$SMQ020==2]<-'Never'; d$smoking[d$SMQ020==1 & d$SMQ040==3]<-'Former'; d$smoking[d$SMQ020==1 & d$SMQ040%in%c(1,2)]<-'Current'; d$smoking<-relevel(factor(d$smoking),'Never')
d$alcohol<-NA_character_; d$alcohol[d$ALQ101==1]<-'Past-year drinker'; d$alcohol[d$ALQ101==2 & d$ALQ110==1]<-'Former/light lifetime'; d$alcohol[d$ALQ101==2 & d$ALQ110==2]<-'Lifetime <12 drinks'; d$alcohol<-relevel(factor(d$alcohol),'Lifetime <12 drinks')
d$phys_active<-ifelse(d$PAQ650==1|d$PAQ665==1,1,ifelse(d$PAQ650==2&d$PAQ665==2,0,NA)); d$diabetes<-ifelse(d$DIQ010==1,1,ifelse(d$DIQ010%in%c(2,3),0,NA)); d$hypertension<-ifelse(d$BPQ020==1,1,ifelse(d$BPQ020==2,0,NA))
hv<-c('MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160F'); d$cvd<-apply(d[,hv],1,function(x)if(any(x==1,na.rm=TRUE))1 else if(all(x==2))0 else NA); d$cancer<-ifelse(d$MCQ220==1,1,ifelse(d$MCQ220==2,0,NA))
f1<-Surv(time,event)~vcf_grade2+low_grip+age10+female+race
f2<-update(f1,.~.+education+INDFMPIR)
f3<-update(f2,.~.+bmi+smoking+alcohol+phys_active)
f4<-update(f3,.~.+diabetes+hypertension+cvd+cancer)
allv<-unique(c(all.vars(f4),'WTMEC2YR','SDMVSTRA','SDMVPSU'))
cc<-d[complete.cases(d[,allv]),]; stopifnot(nrow(cc)==2623,sum(cc$event)==220)
des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC2YR,nest=TRUE,data=cc)
# survey 4.5 attempts an omnibus Wald inversion after fitting. The fully adjusted
# model has more coefficients than design df, so that omnibus matrix is singular.
# Replace only that nonessential global-test assignment; coefficients and the
# package's svyrecvar Taylor-linearized covariance calculation are unchanged.
.svyfun <- get('svycoxph.survey.design', envir=asNamespace('survey'))
.svylines <- deparse(.svyfun)
.svylines[grepl('wald.test', .svylines, fixed=TRUE)] <- '        g$wald.test <- NA_real_'
.svyfun2 <- eval(parse(text=.svylines)[[1]], envir=asNamespace('survey'))
assignInNamespace('svycoxph.survey.design', .svyfun2, ns='survey')
registerS3method('svycoxph','survey.design',.svyfun2,envir=asNamespace('survey'))
extract<-function(fit,label){b<-coef(fit); V<-vcov(fit); se<-sqrt(diag(V)); data.frame(model=label,N=nrow(cc),deaths=sum(cc$event),term=names(b),HR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=2*pnorm(-abs(b/se)),SE=se)}
fits<-list(); for(i in seq_along(list(f1,f2,f3,f4))){cat('fit',i,'\\n'); fits[[i]]<-svycoxph(list(f1,f2,f3,f4)[[i]],des)}; labs<-c('Model 1: demographics','Model 2: + socioeconomic','Model 3: + BMI and lifestyle','Model 4: + comorbidities')
seqall<-do.call(rbind,Map(extract,fits,labs)); seqmain<-seqall[seqall$term%in%c('vcf_grade2','low_grip'),]
# flexible continuous covariates, 4 df for age and 3 df for BMI/PIR
f_age<-Surv(time,event)~vcf_grade2+low_grip+ns(RIDAGEYR,df=4)+female+race+education+INDFMPIR+bmi+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer
f_allflex<-Surv(time,event)~vcf_grade2+low_grip+ns(RIDAGEYR,df=4)+female+race+education+ns(INDFMPIR,df=3)+ns(bmi,df=3)+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer
flex<-rbind(extract(svycoxph(f_age,des),'Flexible age (natural cubic spline, 4 df)'),extract(svycoxph(f_allflex,des),'Flexible age, BMI and PIR'))
flexmain<-flex[flex$term%in%c('vcf_grade2','low_grip'),]
# original custom Taylor implementation on identical data/formula
custom<-function(formula,data,label){z<-data; z$w_norm<-z$WTMEC2YR/mean(z$WTMEC2YR); fit<-coxph(formula,data=z,weights=w_norm,ties='breslow',x=TRUE); db<-residuals(fit,type='dfbeta',weighted=TRUE); agg<-rowsum(db,interaction(z$SDMVSTRA,z$SDMVPSU,drop=TRUE),reorder=FALSE); st<-as.numeric(sub('\\..*$','',rownames(agg))); V<-matrix(0,ncol(db),ncol(db)); for(h in unique(st)){u<-agg[st==h,,drop=FALSE];m<-nrow(u);if(m>1){uc<-sweep(u,2,colMeans(u));V<-V+m/(m-1)*crossprod(uc)}}; b<-coef(fit);se<-sqrt(diag(V));data.frame(method=label,term=names(b),beta=b,SE=se,HR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se))}
cu<-custom(f4,cc,'Custom Taylor'); sv<-extract(svycoxph(f4,des),'survey::svycoxph'); sv2<-data.frame(method='survey::svycoxph',term=sv$term,beta=log(sv$HR),SE=sv$SE,HR=sv$HR,CI_low=sv$CI_low,CI_high=sv$CI_high)
val<-merge(cu,sv2,by='term',suffixes=c('_custom','_survey')); val$beta_abs_diff<-abs(val$beta_custom-val$beta_survey); val$SE_abs_diff<-abs(val$SE_custom-val$SE_survey); val$SE_rel_diff_percent<-100*val$SE_abs_diff/val$SE_survey
# output main exposures and validation
write.csv(seqmain,file.path(out,'four_corrections_common_sample_models.csv'),row.names=FALSE)
write.csv(flexmain,file.path(out,'four_corrections_flexible_covariates.csv'),row.names=FALSE)
write.csv(val,file.path(out,'four_corrections_svycox_validation.csv'),row.names=FALSE)
meta<-data.frame(item=c('R version','survey version','survival version','design degrees of freedom','strata','PSUs','lonely PSU option','N','deaths'),value=c(R.version.string,as.character(packageVersion('survey')),as.character(packageVersion('survival')),degf(des),length(unique(cc$SDMVSTRA)),nrow(unique(cc[c('SDMVSTRA','SDMVPSU')])),getOption('survey.lonely.psu'),nrow(cc),sum(cc$event)))
write.csv(meta,file.path(out,'four_corrections_analysis_metadata.csv'),row.names=FALSE)
# formatted Table 2
fm<-seqmain; fm$Exposure<-ifelse(fm$term=='vcf_grade2','Grade >=2 VCF','FNIH low grip strength'); fm$`HR (95% CI)`<-sprintf('%.2f (%.2f-%.2f)',fm$HR,fm$CI_low,fm$CI_high); fm$`P value`<-ifelse(fm$p<.001,'<0.001',sprintf('%.3f',fm$p)); tab2<-fm[,c('model','N','deaths','Exposure','HR (95% CI)','P value')]; names(tab2)[1:3]<-c('Model','N','Deaths'); write.csv(tab2,file.path(out,'Table2_common_sample_svycoxph.csv'),row.names=FALSE)
# validation/spline sensitivity table
sm<-rbind(data.frame(Analysis=flexmain$model,flexmain[,c('term','N','deaths','HR','CI_low','CI_high','p')]),data.frame(Analysis='Standard survey::svycoxph validation',sv[sv$term%in%c('vcf_grade2','low_grip'),c('term','N','deaths','HR','CI_low','CI_high','p')]))
sm$Exposure<-ifelse(sm$term=='vcf_grade2','Grade >=2 VCF','FNIH low grip strength'); sm$`HR (95% CI)`<-sprintf('%.2f (%.2f-%.2f)',sm$HR,sm$CI_low,sm$CI_high); sm$`P value`<-ifelse(sm$p<.001,'<0.001',sprintf('%.3f',sm$p)); write.csv(sm[,c('Analysis','Exposure','N','deaths','HR','CI_low','CI_high','p','HR (95% CI)','P value')],file.path(out,'Table3_four_corrections_additions.csv'),row.names=FALSE)
# common-sample forest figure
pm<-transform(seqmain,Exposure=ifelse(term=='vcf_grade2','Grade >=2 VCF','FNIH low grip strength'),Model=factor(model,levels=labs)); pd<-position_dodge(width=.55)
p<-ggplot(pm,aes(x=HR,y=Model,color=Exposure,shape=Exposure))+geom_vline(xintercept=1,linetype=2,color='grey55')+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.18,position=pd,linewidth=.65)+geom_point(position=pd,size=2.6)+scale_x_log10(breaks=c(.5,1,2,4),limits=c(.5,5))+scale_color_manual(values=c('Grade >=2 VCF'='#0072B2','FNIH low grip strength'='#D55E00'))+labs(x='Hazard ratio (95% CI)',y=NULL,color=NULL,shape=NULL)+theme_classic(base_size=11)+theme(legend.position='bottom',axis.text.y=element_text(color='black'))
ggsave(file.path(out,'Figure_common_sample_sequential_cox.png'),p,width=7.2,height=4.6,dpi=400,bg='white'); ggsave(file.path(out,'Figure_common_sample_sequential_cox.tiff'),p,width=7.2,height=4.6,dpi=600,compression='lzw',bg='white')
print(seqmain); print(flexmain); print(val[val$term%in%c('vcf_grade2','low_grip'),]); print(meta)




