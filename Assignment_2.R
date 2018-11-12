####### Assignment 2 #############
library(survival)
library(foreign)
library(splines)

d=read.spss('renaltx.sav',to.data.frame=TRUE)
apply(is.na(d),2,sum)
min(d$gsurv)
#categorical: aantalre uprotein cregsh, binary: dgf, continuous:creat prac predias acclft

#check linearity and proportionality to see if you can do cox
#transform data?
#otherwise do splines
#polynomial?

#Kaplan Meier for survival months and survival status
plot(survfit(Surv(gsurv,gstatus)~1,data=d),xlab = 'Months since start', ylab = 'Survival', main = 'Kaplan-Meier')

#logrank for two groups (dgf)
plot(survfit(Surv(gsurv,gstatus)~dgf,data=d),col=c(2,3))
survdiff(Surv(gsurv,gstatus)~dgf,data=d)

#Univariate coxgression for all covariates
univar1=coxph(Surv(gsurv,gstatus)~aantalre,data=d)
plot(cox.zph(univar1))
summary(univar1)
#als factor als je niet blij bent met linearity
univar1a=coxph(Surv(gsurv,gstatus)~as.factor(aantalre),data=d)
plot(cox.zph(univar1a))
summary(univar1a)
#linearity kan gecheckt worden bij coefficients
#de niet lineare factor (not as.factor) geeft coef = 2, dus dan verwacht je level 1 is 2, level 2 is 4 level 3 is 8 etc
#dit kan je dan checken met de coeff voor het as.factor model


#dan proportionality checken door cox.zph te plotten, is de concordance uit het summary van het model goed
#geschat als baseline voor de proportionality (we verwachten een rechte lijn) =0.7 voor niet as.factor
#concordance lijkt hetzelfde voor beide modellen
#wat te doen als niet proportional??? google


#Linearity check voor de continuous variables
#als niet linear dan splines om dat te fixen
#c statistic tussen modellen

univar2=coxph(Surv(gsurv,gstatus)~creat,data=d)
summary(univar2)
plot(cox.zph(univar2), main='proportionality of creatinine', xlab = 'n.o. months',ylab = 'corresponding exp of singular coefficient')
abline(h=0)



univar2a=coxph(Surv(gsurv,gstatus)~pspline(creat, df=10, caic=T),data=d)
summary(univar2a)
#linear model maken en kijken of verworpen wordt
hh2 = predict(univar2, type = 'terms',se.fit=T)
hh2a=predict(univar2a, type="terms", se.fit=T)
plot(d$creat,hh2$fit)
plot(d$creat,hh2a$fit,main = 'linearity of creatinine',xlab = 'actual creatinine conc.', ylab = 'predicted creat conc.')


#univariate list of models for all covariates
col_names=colnames(d)
N <- length(col_names)
univarlist <- vector("list", N)
for(i in 1:N) {
  indexformula <- as.formula(sprintf("Surv(gsurv, gstatus) ~ %s", col_names[i]))
  univar <-  coxph(indexformula,data=d,na.action=na.omit)
  univarlist[[i]] <- univar
}
setNames(univarlist, paste("univar", 1:length(col_names),sep=''))
univarlist[[1]]

AIC(univarlist[[10]])

#check linearity assumption of quantative
univar1a=coxph(Surv(gsurv,gstatus)~pspline(acclft,df=0),data=d)
univar2a=coxph(Surv(gsurv,gstatus)~pspline(predias,df=0),data=d)
univar3a=coxph(Surv(gsurv,gstatus)~pspline(creat,df=0),data=d)
univar4a=coxph(Surv(gsurv,gstatus)~pspline(prac,df=0),data=d)

hh1a=predict(univar1a, type="terms", se.fit=T)
plot(d$acclft[is.na(d$acclft)==FALSE],hh1a$fit)

hh2a=predict(univar2a, type="terms", se.fit=T)
plot(d$predias[is.na(d$predias)==FALSE],hh2a$fit)

hh3a=predict(univar3a, type="terms", se.fit=T)
plot(d$creat[is.na(d$creat)==FALSE],hh3a$fit)
#Plot creat against survival time
plot(d$creat,d$gstatus)

hh4a=predict(univar4a, type="terms", se.fit=T)
plot(d$prac[is.na(d$prac)==FALSE],hh4a$fit)

#Find best splines for quantative data, variable aantalre
AIClist = data.frame(rep(0,15))
for (i in 1:15){
  AIClist[i,] <- AIC(coxph(Surv(gsurv,gstatus)~ns(acclft,df=i),data=d))
}
plot(1:length(AIClist$rep.0..15.),AIClist$rep.0..15., main='AIC versus d.o.f.',xlab = 'd.o.f.',ylab = 'd.o.f.')

#model comparisons
univar5a=coxph(Surv(gsurv,gstatus)~prac + pspline(creat,df=10) + predias + acclft + dgf + aantalre + uprotein + cregsh,data=d)
univar6a=coxph(Surv(gsurv,gstatus)~pspline(creat,df=10)+ acclft + uprotein + dgf,data=d)
anova(univar5a)
anova(univar6a)
AIC(univar5a)
AIC(univar6a)

#binning and factorising prac
pracbin <- cut(d$prac, c(0,20, 40,80, 120, max(d$prac)), include.lowest = TRUE, labels=c("low", "med", "high","higher","higherst"))
d$pracbin <- pracbin
univar5a=coxph(Surv(gsurv,gstatus)~as.factor(pracbin),data=d)
predict(univar5a, type="terms", se.fit=T)
plot(d$prac[is.na(d$pracbin)==FALSE],hh5a$fit)
