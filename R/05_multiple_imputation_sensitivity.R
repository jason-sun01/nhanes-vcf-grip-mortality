library(survival)
library(nnet)

set.seed(20260901)
d0 <- read.csv("data_processed/nhanes_full_covariate_analysis_dataset.csv",na.strings=c("NA",""))
d0$race <- factor(d0$RIDRETH1)
d0$education <- factor(d0$education,levels=1:5)
d0$smoking <- factor(d0$smoking,levels=c("Never","Former","Current"))
d0$alcohol <- factor(d0$alcohol,levels=c("Lifetime <12 drinks","Former/light lifetime","Past-year drinker"))
d0$logtime <- log(pmax(d0$time,1/12))

impvars <- c("education","INDFMPIR","bmi","smoking","alcohol","phys_active",
             "diabetes","hypertension","cvd","cancer")
missing0 <- lapply(d0[impvars],is.na); names(missing0)<-impvars

initial_fill <- function(x){
  miss<-is.na(x); obs<-x[!miss]
  if(is.factor(x)) x[miss]<-sample(obs,sum(miss),replace=TRUE)
  else x[miss]<-sample(obs,sum(miss),replace=TRUE)
  x
}

pmm_update <- function(dat,var,formula,k=5){
  mis<-missing0[[var]]; if(!any(mis)) return(dat[[var]])
  obs<-!mis; boot<-sample(which(obs),sum(obs),replace=TRUE)
  fit<-lm(formula,data=dat[boot,])
  po<-predict(fit,newdata=dat[obs,]); pm<-predict(fit,newdata=dat[mis,])
  donors<-which(obs); out<-dat[[var]]
  for(j in seq_along(pm)){
    near<-order(abs(po-pm[j]))[seq_len(min(k,length(po)))]
    out[which(mis)[j]]<-dat[[var]][sample(donors[near],1)]
  }
  out
}

cat_update <- function(dat,var,formula){
  mis<-missing0[[var]]; if(!any(mis)) return(dat[[var]])
  obs<-!mis; boot<-sample(which(obs),sum(obs),replace=TRUE)
  fit<-multinom(formula,data=dat[boot,],trace=FALSE,MaxNWts=5000)
  pr<-predict(fit,newdata=dat[mis,],type="probs")
  lev<-levels(dat[[var]])
  if(is.null(dim(pr))){
    if(sum(mis)==1 && length(pr)==length(lev)) pr<-matrix(pr,nrow=1,dimnames=list(NULL,lev))
    else pr<-cbind(1-pr,pr)
  }
  if(is.null(colnames(pr))) colnames(pr)<-lev
  out<-dat[[var]]; idx<-which(mis)
  for(j in seq_along(idx)){
    pp<-pr[j,lev]; pp[is.na(pp)]<-0; pp<-pp/sum(pp)
    out[idx[j]]<-sample(lev,1,prob=pp)
  }
  factor(out,levels=lev)
}

bin_update <- function(dat,var,formula){
  mis<-missing0[[var]]; if(!any(mis)) return(dat[[var]])
  obs<-!mis; boot<-sample(which(obs),sum(obs),replace=TRUE)
  fit<-glm(formula,data=dat[boot,],family=binomial())
  pr<-pmin(pmax(predict(fit,newdata=dat[mis,],type="response"),.001),.999)
  out<-dat[[var]]; out[mis]<-rbinom(sum(mis),1,pr); out
}

cox_design_two <- function(dat){
  dat$w_norm<-dat$WTMEC2YR/mean(dat$WTMEC2YR)
  fit<-coxph(Surv(time,event)~vcf_grade2+low_grip+age10+female+race+education+
      INDFMPIR+bmi+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer,
      data=dat,weights=w_norm,ties="breslow",x=TRUE)
  db<-residuals(fit,type="dfbeta",weighted=TRUE)
  agg<-rowsum(db,interaction(dat$SDMVSTRA,dat$SDMVPSU,drop=TRUE),reorder=FALSE)
  stratum<-as.numeric(sub("\\..*$","",rownames(agg))); V<-matrix(0,ncol(db),ncol(db))
  for(h in unique(stratum)){u<-agg[stratum==h,,drop=FALSE];m<-nrow(u);if(m>1){uc<-sweep(u,2,colMeans(u));V<-V+m/(m-1)*crossprod(uc)}}
  ix<-match(c("vcf_grade2","low_grip"),names(coef(fit)))
  data.frame(term=c("vcf_grade2","low_grip"),beta=coef(fit)[ix],variance=diag(V)[ix])
}

M<-20; fits<-vector("list",M)
for(m in seq_len(M)){
  z<-d0; for(v in impvars) z[[v]]<-initial_fill(z[[v]])
  for(cycle in 1:3){
    z$INDFMPIR<-pmm_update(z,"INDFMPIR",INDFMPIR~age10+female+race+education+bmi+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer+vcf_grade2+low_grip+event+logtime)
    z$bmi<-pmm_update(z,"bmi",bmi~age10+female+race+education+INDFMPIR+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer+vcf_grade2+low_grip+event+logtime)
    z$education<-cat_update(z,"education",education~age10+female+race+INDFMPIR+bmi+smoking+alcohol+event+logtime)
    z$smoking<-cat_update(z,"smoking",smoking~age10+female+race+education+INDFMPIR+bmi+alcohol+event+logtime)
    z$alcohol<-cat_update(z,"alcohol",alcohol~age10+female+race+education+INDFMPIR+bmi+smoking+event+logtime)
    z$phys_active<-bin_update(z,"phys_active",phys_active~age10+female+race+education+INDFMPIR+bmi+event+logtime)
    z$diabetes<-bin_update(z,"diabetes",diabetes~age10+female+race+bmi+event+logtime)
    z$hypertension<-bin_update(z,"hypertension",hypertension~age10+female+race+bmi+event+logtime)
    z$cvd<-bin_update(z,"cvd",cvd~age10+female+race+bmi+event+logtime)
    z$cancer<-bin_update(z,"cancer",cancer~age10+female+race+bmi+event+logtime)
  }
  fits[[m]]<-cox_design_two(z); fits[[m]]$imputation<-m
}
allfit<-do.call(rbind,fits); write.csv(allfit,"results/nhanes_MI_20_imputation_estimates.csv",row.names=FALSE)

pool_one<-function(x){
  qbar<-mean(x$beta); ubar<-mean(x$variance); B<-var(x$beta); T<-ubar+(1+1/M)*B
  r<-((1+1/M)*B)/ubar; df<-(M-1)*(1+1/r)^2
  se<-sqrt(T); p<-2*pt(-abs(qbar/se),df=df)
  data.frame(term=x$term[1],M=M,logHR=qbar,SE=se,df=df,HR=exp(qbar),
             CI_low=exp(qbar-qt(.975,df)*se),CI_high=exp(qbar+qt(.975,df)*se),p=p,
             fraction_missing_information=(r+2/(df+3))/(r+1))
}
pooled<-do.call(rbind,lapply(split(allfit,allfit$term),pool_one))
write.csv(pooled,"results/nhanes_MI_20_pooled_results.csv",row.names=FALSE)
print(pooled)

