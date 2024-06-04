# Weighted code


hcctcga <- read.csv2("livercancerTCGAclinical.csv",header=T,sep=",")
hccbclc <- read.csv("HCCBCLC.csv",header=T)
names(hccbclc)[8] <- "age"
names(hccbclc)[9] <- "gender"
hccbclc$gender <- ifelse(hccbclc$gender == 1,"FEMALE","MALE")
hcctcga$gender <- as.factor(hcctcga$gender)
hccbclc$gender <- as.factor(hccbclc$gender)
hccbclc$futime <- hccbclc$OS
hccbclc$fustat <- hccbclc$D5
library(survival)
library(randomForestSRC)
library(survminer)

fit <- survfit(Surv(futime,fustat)~gender,data=hcctcga)
ggsurvplot(fit,risk.table=T,tables.height=0.3)

fit <- survfit(Surv(futime,fustat)~gender,data=hccbclc)
ggsurvplot(fit,risk.table=T,tables.height=0.3)



cph1 <- coxph(Surv(futime,fustat)~age+gender,data=hcctcga)
wts <- predict(cph1,newdata=hccbclc,type="risk")

cph2 <- coxph(Surv(OS,D5)~BCLC+ER+LR,data=hccbclc,weights = wts)
# Now with no weights
cph3 <- coxph(Surv(OS,D5)~BCLC+ER+LR,data=hccbclc)
cph4 <- coxph(Surv(OS,D5)~BCLC+ER+LR+tmp,data=hccbclc)
n <- length(wts)
# Try new approach
B <- 1000
coefs <- matrix(0,nrow=B,ncol=3)
for (i in 1:B) {
  set.seed(520+i)
  simy <- rexp(n)
  pwts <- simy*wts
  tmpcph <-  coxph(Surv(OS,D5)~BCLC+ER+LR,data=hccbclc,weights = pwts)
  coefs[i,] <- coef(tmpcph)
  cat(i,"\n")
}


# Now try random forests

obj <- rfsrc(Surv(futime,fustat)~age+gender,data=hcctcga,ntree=1000)
wts <- predict(obj,hccbclc)$chf[,120]

cph2 <- coxph(Surv(OS,D5)~BCLC+ER+LR,data=hccbclc,weights = wts)
# Try new approach
B <- 1000
coefs <- matrix(0,nrow=B,ncol=3)
n <- length(wts)
for (i in 1:B) {
  set.seed(520+i)
  simy <- rexp(n)
  pwts <- simy*wts
  tmpcph <-  coxph(Surv(OS,D5)~BCLC+ER+LR,data=hccbclc,weights = pwts)
  coefs[i,] <- coef(tmpcph)
  cat(i,"\n")
}

# next, additive hazards with proportional hazards weights
ah1 <- aalen(Surv(OS,D5)~const(BCLC)+const(ER)+const(LR),data=hccbclc)

cph1 <- coxph(Surv(futime,fustat)~age+gender,data=hcctcga)
wts <- predict(cph1,newdata=hccbclc,type="risk")

ah2 <- aalen(Surv(OS,D5)~const(BCLC)+const(ER)+const(LR),data=hccbclc,weights=wts)
# Try new approach
B <- 1000
coefs <- matrix(0,nrow=B,ncol=3)
n <- length(wts)
for (i in 1:B) {
  set.seed(520+i)
  simy <- rexp(n)
  pwts <- simy*wts
  tmpah<-  aalen(Surv(OS,D5)~const(BCLC)+const(ER)+const(LR),
                 data=hccbclc,weights=pwts)
  coefs[i,] <- coef(tmpah)[,1]
  cat(i,"\n")
}

# Next, additive hazards with random forest weights
obj <- rfsrc(Surv(futime,fustat)~age+gender,data=hcctcga,ntree=1000)
wts <- predict(obj,hccbclc)$chf[,120]
ah3 <- aalen(Surv(OS,D5)~const(BCLC)+const(ER)+const(LR),data=hccbclc,weights=wts)
# Try new approach
B <- 1000
coefs <- matrix(0,nrow=B,ncol=3)
n <- length(wts)
for (i in 1:B) {
  set.seed(520+i)
  simy <- rexp(n)
  pwts <- simy*wts
  tmpah<-  aalen(Surv(OS,D5)~const(BCLC)+const(ER)+const(LR),
                 data=hccbclc,weights=pwts)
  coefs[i,] <- coef(tmpah)[,1]
  cat(i,"\n")
}


