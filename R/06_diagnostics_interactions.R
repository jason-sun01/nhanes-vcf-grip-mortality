library(survival)

d <- read.csv("data_processed/nhanes_full_covariate_analysis_dataset.csv",na.strings=c("NA",""))
d$race <- factor(d$RIDRETH1)
d$education <- factor(d$education,levels=1:5)
d$smoking <- factor(d$smoking,levels=c("Never","Former","Current"))
d$alcohol <- factor(d$alcohol,levels=c("Lifetime <12 drinks","Former/light lifetime","Past-year drinker"))
d$any_vcf <- as.integer(d$DXXVFAST==2)
d$grip_lower5 <- -d$grip_max/5
d$older65 <- as.integer(d$RIDAGEYR>=65)

base_cov <- "+age10+female+race+education+INDFMPIR+bmi+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer"
f4 <- as.formula(paste("Surv(time,event)~vcf_grade2+low_grip",base_cov))
vars <- all.vars(f4)
cc <- d[complete.cases(d[,unique(c(vars,"WTMEC2YR","SDMVSTRA","SDMVPSU"))]),]

design_fit <- function(formula,data,label){
  z <- data[complete.cases(data[,unique(c(all.vars(formula),"WTMEC2YR","SDMVSTRA","SDMVPSU"))]),]
  z$w_norm <- z$WTMEC2YR/mean(z$WTMEC2YR)
  fit <- coxph(formula,data=z,weights=w_norm,ties="breslow",x=TRUE)
  db <- residuals(fit,type="dfbeta",weighted=TRUE)
  agg <- rowsum(db,interaction(z$SDMVSTRA,z$SDMVPSU,drop=TRUE),reorder=FALSE)
  stratum <- as.numeric(sub("\\..*$","",rownames(agg))); V <- matrix(0,ncol(db),ncol(db))
  for(h in unique(stratum)){u<-agg[stratum==h,,drop=FALSE];m<-nrow(u);if(m>1){uc<-sweep(u,2,colMeans(u));V<-V+m/(m-1)*crossprod(uc)}}
  b<-coef(fit);se<-sqrt(diag(V))
  tab<-data.frame(analysis=label,N=nrow(z),deaths=sum(z$event),term=names(b),HR=exp(b),
    CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=2*pnorm(-abs(b/se)),row.names=NULL)
  list(table=tab,fit=fit,data=z,dfbeta=db)
}

primary <- design_fit(f4,cc,"Primary complete-case Model 4")

# Proportional-hazards diagnostic. cox.zph assesses coefficient-time trends;
# survey-weighted model fit is retained, but these P values are diagnostic.
zph <- cox.zph(primary$fit,transform="km")
zph_tab <- data.frame(term=rownames(zph$table),chisq=zph$table[,"chisq"],
                      df=zph$table[,"df"],p=zph$table[,"p"],row.names=NULL)
write.csv(zph_tab,"results/nhanes_PH_assumption_tests.csv",row.names=FALSE)

# Landmark analysis: condition on surviving beyond two years and reset time origin.
land <- cc[cc$time>2,]; land$time <- land$time-2
landmark <- design_fit(f4,land,"Two-year landmark")

# Alternative definitions and continuous grip.
f_any <- as.formula(paste("Surv(time,event)~any_vcf+low_grip",base_cov))
f_cont <- as.formula(paste("Surv(time,event)~vcf_grade2+grip_lower5",base_cov))
anygrade <- design_fit(f_any,cc,"Any-grade VCF")
continuous <- design_fit(f_cont,cc,"Grip per 5-kg lower")

# Interaction screens. Main effects stay in each model.
f_sex <- as.formula(paste("Surv(time,event)~vcf_grade2*female+low_grip*female",sub("+female","",base_cov,fixed=TRUE)))
f_age <- as.formula(paste("Surv(time,event)~vcf_grade2*older65+low_grip*older65",base_cov))
f_joint <- as.formula(paste("Surv(time,event)~vcf_grade2*low_grip",base_cov))
sexint <- design_fit(f_sex,cc,"Sex interactions")
ageint <- design_fit(f_age,cc,"Age >=65 interactions")
jointint <- design_fit(f_joint,cc,"VCF x low-grip interaction")

# Sex-stratified estimates to interpret the statistically notable VCF-by-sex term.
cov_nosex <- sub("+female","",base_cov,fixed=TRUE)
f_strat <- as.formula(paste("Surv(time,event)~vcf_grade2+low_grip",cov_nosex))
menfit <- design_fit(f_strat,cc[cc$female==0,],"Men")
womenfit <- design_fit(f_strat,cc[cc$female==1,],"Women")
sex_counts <- aggregate(event~female+vcf_grade2,cc,function(x)c(N=length(x),deaths=sum(x)))
sex_counts <- data.frame(female=sex_counts$female,vcf_grade2=sex_counts$vcf_grade2,
                         N=sex_counts$event[,"N"],deaths=sex_counts$event[,"deaths"])
write.csv(sex_counts,"results/nhanes_sex_vcf_event_counts.csv",row.names=FALSE)

# Influence audit and deletion sensitivity using union of top five absolute
# weighted dfbeta values for each primary exposure.
db <- primary$dfbeta; colnames(db)<-names(coef(primary$fit))
top <- unique(c(order(abs(db[,"vcf_grade2"]),decreasing=TRUE)[1:5],
                order(abs(db[,"low_grip"]),decreasing=TRUE)[1:5]))
influence <- data.frame(SEQN=primary$data$SEQN[top],event=primary$data$event[top],
  vcf_grade2=primary$data$vcf_grade2[top],low_grip=primary$data$low_grip[top],
  dfbeta_vcf=db[top,"vcf_grade2"],dfbeta_low_grip=db[top,"low_grip"])
write.csv(influence,"results/nhanes_top_influence_observations.csv",row.names=FALSE)
deleted <- design_fit(f4,primary$data[-top,],"Exclude top influence observations")

alltab <- rbind(primary$table,landmark$table,anygrade$table,continuous$table,
                sexint$table,ageint$table,jointint$table,menfit$table,womenfit$table,deleted$table)
keepterms <- c("vcf_grade2","low_grip","any_vcf","grip_lower5",
               "vcf_grade2:female","female:low_grip",
               "vcf_grade2:older65","older65:low_grip","vcf_grade2:low_grip")
key <- alltab[alltab$term %in% keepterms,]
write.csv(alltab,"results/nhanes_diagnostics_sensitivity_all_terms.csv",row.names=FALSE)
write.csv(key,"results/nhanes_diagnostics_sensitivity_key_results.csv",row.names=FALSE)
print(zph_tab); print(key)

