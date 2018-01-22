library(ggplot2)
library(RankProd)
library(samr)
library(reshape2)
library(gplots)
library(coin)
library(metap)
library(igraph)
library(plotly)



feng=read.csv("feng.csv")
zeller=read.csv("zeller.csv")
yu=read.csv("yu.csv")
vogtmann=read.csv("vogtmann.csv")
taxonomy=read.csv("Species Taxonomy.csv")


feng=feng[feng$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(feng$age,feng$BMI,feng$gender)
feng=feng[complete.cases(phenotype),]
feng$grouped_diagnosis=as.character(feng$grouped_diagnosis)
feng.abundance=feng[,-c(1:14)]
feng.abundance=feng.abundance[,names(feng.abundance) %in% taxonomy$OTUcode]
x=(feng.abundance>0)*1
sum=colSums(x)
feng.abundance=feng.abundance[,(sum>0.5*nrow(feng.abundance))]



zeller=zeller[zeller$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(zeller$age,zeller$BMI,zeller$gender)
zeller=zeller[complete.cases(phenotype),]
zeller$grouped_diagnosis=as.character(zeller$grouped_diagnosis)
zeller.abundance=zeller[,-c(1:14)]
zeller.abundance=zeller.abundance[,names(zeller.abundance) %in% taxonomy$OTUcode]
x=(zeller.abundance>0)*1
sum=colSums(x)
zeller.abundance=zeller.abundance[,(sum>0.5*nrow(zeller.abundance))]




yu=yu[yu$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(yu$age,yu$BMI,yu$gender)
yu=yu[complete.cases(phenotype),]
yu$grouped_diagnosis=as.character(yu$grouped_diagnosis)
yu.abundance=yu[,-c(1:14)]
yu.abundance=yu.abundance[,names(yu.abundance) %in% taxonomy$OTUcode]
x=(yu.abundance>0)*1
sum=colSums(x)
yu.abundance=yu.abundance[,(sum>0.5*nrow(yu.abundance))]



vogtmann=vogtmann[vogtmann$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(vogtmann$age,vogtmann$BMI,vogtmann$gender)
vogtmann=vogtmann[complete.cases(phenotype),]
vogtmann$grouped_diagnosis=as.character(vogtmann$grouped_diagnosis)
vogtmann.abundance=vogtmann[,-c(1:14)]
vogtmann.abundance=vogtmann.abundance[,names(vogtmann.abundance) %in% taxonomy$OTUcode]
x=(vogtmann.abundance>0)*1
sum=colSums(x)
vogtmann.abundance=vogtmann.abundance[,(sum>0.5*nrow(vogtmann.abundance))]




######################### Select Differentially Abundant Bacteria ################################

### Adjust for confounders

for (i in 1:ncol(yu.abundance)){
  y=log2(yu.abundance[,i]+1)
  lm.fit=lm(y~age,data = yu)
  s=summary(lm.fit)
  p.value=s$coefficients[2,4]
  if (p.value<0.05){
    yu.abundance[,i]=2^(lm.fit$residuals)
  }
}

for (i in 1:ncol(zeller.abundance)){
  y=log2(zeller.abundance[,i]+1)
  lm.fit=lm(y~age,data = zeller)
  s=summary(lm.fit)
  p.value=s$coefficients[2,4]
  if (p.value<0.05){
    zeller.abundance[,i]=2^(lm.fit$residuals)
  }
}

### Select Bacteria with same change direction
feng.1=feng.abundance[feng$grouped_diagnosis=="normal_control",]
feng.2=feng.abundance[feng$grouped_diagnosis=="CRC",]
zeller.1=zeller.abundance[zeller$grouped_diagnosis=="normal_control",]
zeller.2=zeller.abundance[zeller$grouped_diagnosis=="CRC",]
yu.1=yu.abundance[yu$grouped_diagnosis=="normal_control",]
yu.2=yu.abundance[yu$grouped_diagnosis=="CRC",]
vogtmann.1=vogtmann.abundance[vogtmann$grouped_diagnosis=="normal_control",]
vogtmann.2=vogtmann.abundance[vogtmann$grouped_diagnosis=="CRC",]


ch.1=apply(feng.2,2,median)/apply(feng.1,2,median)
ch.1=data.frame(ch.1,OTU=names(feng.abundance))
ch.2=apply(zeller.2,2,median)/apply(zeller.1,2,median)
ch.2=data.frame(ch.2,OTU=names(zeller.abundance))
ch.3=apply(yu.2,2,median)/apply(yu.1,2,median)
ch.3=data.frame(ch.3,OTU=names(yu.abundance))
ch.4=apply(vogtmann.2,2,median)/apply(vogtmann.1,2,median)
ch.4=data.frame(ch.4,OTU=names(vogtmann.abundance))

change=merge(ch.1,ch.2,by="OTU")
change=merge(change,ch.3,by="OTU")
change=merge(change,ch.4,by="OTU")
change=change[complete.cases(change),]

### Test the selected bacteria with RankSum 

feng.abundance=t(feng.abundance)
names(feng.abundance)=feng$Sample.ID
zeller.abundance=t(zeller.abundance)
names(zeller.abundance)=zeller$Sample.ID
yu.abundance=t(yu.abundance)
names(yu.abundance)=yu$Sample.ID
vogtmann.abundance=t(vogtmann.abundance)
names(vogtmann.abundance)=vogtmann$Sample.ID



feng.abundance=feng.abundance[(row.names(feng.abundance) %in% change$OTU),]
zeller.abundance=zeller.abundance[(row.names(zeller.abundance) %in% change$OTU),]
yu.abundance=yu.abundance[(row.names(yu.abundance) %in% change$OTU),]
vogtmann.abundance=vogtmann.abundance[(row.names(vogtmann.abundance) %in% change$OTU),]


abundance=data.frame(feng.abundance,zeller.abundance,yu.abundance,vogtmann.abundance)

origin=c(rep(1,ncol(feng.abundance)),rep(2,ncol(zeller.abundance)),rep(3,ncol(yu.abundance)),rep(4,ncol(vogtmann.abundance)))

x=change[,-1]
x=(x>=1)*1
sum=rowSums(x)

for (i in 1:length(sum)){
  j=0
  if (sum[i]==1){
    j=which(x[i,]==1)
  }
  if (sum[i]==3){
    j=which(x[i,]==0)
  }
  if (j!=0){
    if (j==1){
      m1=wilcox.test(feng.1[,i],feng.2[,i])
      if (m1$p.value<0.05) sum[i]=2
    }
    if (j==2){
      m2=wilcox.test(zeller.1[,i],zeller.2[,i])
      if (m2$p.value<0.05) sum[i]=2
    }
    if (j==3){
      m2=wilcox.test(yu.1[,i],yu.2[,i])
      if (m2$p.value<0.05) sum[i]=2
    }
    if (j==4){
      m2=wilcox.test(vogtmann.1[,i],vogtmann.2[,i])
      if (m2$p.value<0.05) sum[i]=2
    }
  }
}

status1=(feng$grouped_diagnosis=="CRC")*1
status2=(zeller$grouped_diagnosis=="CRC")*1
status3=(yu$grouped_diagnosis=="CRC")*1
status4=(vogtmann$grouped_diagnosis=="CRC")*1
status=c(status1,status2,status3,status4)

selected.abundance=abundance[sum %in% c(3,4),]

output1=RSadvance(selected.abundance,cl=status,origin=origin,num.perm = 100,logged=F)

top.enriched=topGene(output,cutoff=0.01,method="pfp",logged=F,gene.names=row.names(selected.abundance))

selected.abundance=abundance[sum %in% c(0,1),]

output2=RSadvance(selected.abundance,cl=status,origin=origin,num.perm = 100,logged=F)

top.depleted=topGene(output,cutoff=0.01,method="pfp",logged=F,gene.names=row.names(selected.abundance))








#################  Build the Prediction Model With CRC Enriched Species ############################

library(pROC)
library(e1071)
library(plotROC)

marker=c("Otu0302","Otu0918","Otu0309","Otu0291","Otu0915","Otu1887","Otu0314")
feng.marker=feng.abundance[,names(feng.abundance) %in% marker]
zeller.marker=zeller.abundance[,names(zeller.abundance) %in% marker]
yu.marker=yu.abundance[,names(yu.abundance) %in% marker]
vogtmann.marker=vogtmann.abundance[,names(vogtmann.abundance) %in% marker]
cohort=c(rep("feng",nrow(feng.marker)),rep("zeller",nrow(zeller.marker)),
         rep("yu",nrow(yu.marker)),rep("vogtmann",nrow(vogtmann.marker)))

age=c(feng$age,zeller$age,yu$age,vogtmann$age)
gender=c(feng$gender,zeller$gender,yu$gender,vogtmann$gender)
BMI=c(feng$BMI,zeller$BMI,yu$BMI,vogtmann$BMI)
feng.marker=feng.marker/1596464
zeller.marker=zeller.marker/856204
yu.marker=yu.marker/1222507
vogtmann.marker=vogtmann.marker/2419973
marker=rbind(feng.marker,zeller.marker,yu.marker,vogtmann.marker)
marker2=sqrt(marker)
grouped_diagnosis=c(feng$grouped_diagnosis,zeller$grouped_diagnosis,yu$grouped_diagnosis,vogtmann$grouped_diagnosis)

dat_all=data.frame(age,gender,BMI,marker2)


### Tune the parameter through cross-validating AUC

library(caret)
folds=createFolds(c(1:526),k=10)
gamma_list=seq(0.05, 0.5, 0.05)
coef0_list=seq(0, 1.9, 0.1)
degree_list=seq(1, 5, 1)

mean.auc=rep(0,10)
for (i in 1:10){
  for (k in 1:10){
    train.y=grouped_diagnosis[-(folds[[k]])]
    train.x=marker2[-(folds[[k]]),]
    test.y=grouped_diagnosis[(folds[[k]])]
    test.x=marker2[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = gamma_list[i], probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])$auc
    mean.auc[i]=mean.auc[i]+r5/10
  }
}


### Overall Prediction Model with Radial Kernel,gamma=0.2

m1=svm(x=as.matrix(marker2),y=as.factor(grouped_diagnosis), gamma = 0.2, probability = T)
class=predict(m1,marker2,probability = T)
x=attributes(class)
r1=roc(grouped_diagnosis~x$probabilities[,1])

m2=svm(x=as.matrix(dat_all),y=as.factor(grouped_diagnosis), gamma = 0.2, probability = T)
class=predict(m2,dat_all,probability = T)
x=attributes(class)
r2=roc(grouped_diagnosis~x$probabilities[,1])



phenotype=data.frame(age,gender,BMI)

m3=svm(x=as.matrix(phenotype),y=as.factor(grouped_diagnosis), gamma = 0.2, probability = T)
class=predict(m3,phenotype,probability = T)
x=attributes(class)
r3=roc(grouped_diagnosis~x$probabilities[,1])


sensitivity=c(r1$sensitivities,r2$sensitivities,r3$sensitivities)
specificity=c(r1$specificities,r2$specificities,r3$specificities)

cohort=c(rep("Bacteria = 0.80",length(r1$sensitivities)),rep("Clinical = 0.68",length(r2$sensitivities)),
         rep("Bac+Cli = 0.88",length(r3$sensitivities)))
cohort=factor(cohort,levels = c("Bacteria = 0.80","Clinical = 0.68","Bac+Cli = 0.88"))

dat=data.frame(sensitivity,specificity=1-specificity,cohort)

ggplot(aes(x=specificity,y=sensitivity,color=cohort),data=dat)+
  geom_path(size=0.9)+theme_bw()+theme(legend.text = element_text(size=14))+
  scale_color_manual(values=c("#FFCC00","#33CC33","#669999"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+
  theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.65,0.15))+
  theme(legend.title = element_blank())+xlab("1-Specificity")+ylab("Sensitivity")


### Plot the Cross-Validated ROC Curve

roc=rep(list(),400)

mean.auc=0
for (i in 1:40){
  folds=createFolds(c(1:526),k=10)
  for (k in 1:10){
    train.y=grouped_diagnosis[-(folds[[k]])]
    train.x=marker2[-(folds[[k]]),]
    test.y=grouped_diagnosis[(folds[[k]])]
    test.x=marker2[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = 0.2, probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])
    roc[[(i*10+k-10)]]=r5
  }
}

Check_density=function(v,x,y){
  x1=x[1:(length(x)-1)]
  x2=x[2:length(x)]
  return(max(y[(x1<=v)&(x2>=v)]))
}

simple_auc <- function(TPR, FPR){
  # inputs already sorted, worst scores first 
  dFPR <- -c(diff(FPR), 0)
  dTPR <- -c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

mean.auc=0
all.auc=rep(0,length(roc))
for (i in 1:length(roc)){
  mean.auc=mean.auc+roc[[i]]$auc/length(roc)
  all.auc[i]=roc[[i]]$auc
}

roc=roc[order(all.auc)]

specificity=0
for (i in 1:length(roc)){
  specificity=c(specificity,roc[[i]]$specificities)
}
specificity=sort(unique(specificity))
sensitivity=upper=lower=rep(0,length(specificity))
for (i in 1:length(specificity)){
  m=rep(0,length(roc))
  for (k in 1:length(roc)){
    m[k]=Check_density(specificity[i],roc[[k]]$specificities,roc[[k]]$sensitivities)
  }
  sensitivity[i]=mean(m)
  lower[i]=quantile(m,0.025)
  upper[i]=quantile(m,0.975)
}

dat=data.frame(sensitivity,upper,lower,specificity=1-specificity)
dat=rbind(dat,c(0,min(upper),min(lower),0))

ggplot(aes(x=specificity,y=sensitivity),data=dat)+
  geom_path(size=1.1, color="#206999")+theme_bw()+theme(legend.text = element_text(size=14))+
  geom_ribbon(data=dat,aes(x=specificity,ymin=lower,ymax=upper),alpha=0.2)+
  geom_abline(slope=1, intercept=0, linetype=2)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17,face="italic"))+
  theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.7,0.15))+
  theme(legend.title = element_blank())+xlab("1-Specificity")+ylab("Sensitivity")




### Leave One-cohort Out

train.y1=grouped_diagnosis[-c(1:109)]
train.y2=grouped_diagnosis[-c(110:261)]
train.y3=grouped_diagnosis[-c(262:426)]
train.y4=grouped_diagnosis[-c(427:526)]
train.x1=marker2[-c(1:109),]
train.x2=marker2[-c(110:261),]
train.x3=marker2[-c(262:426),]
train.x4=marker2[-c(427:526),]
test.y1=grouped_diagnosis[c(1:109)]
test.y2=grouped_diagnosis[c(110:261)]
test.y3=grouped_diagnosis[c(262:426)]
test.y4=grouped_diagnosis[c(427:526)]
test.x1=marker2[c(1:109),]
test.x2=marker2[c(110:261),]
test.x3=marker2[c(262:426),]
test.x4=marker2[c(427:526),]


folds=createFolds(c(1:417), k=10)
gamma_list=seq(0.1, 1, 0.1)

mean.auc=rep(0,10)
for (i in 1:10){
  for (k in 1:10){
    train.y=train.y1[-(folds[[k]])]
    train.x=train.x1[-(folds[[k]]),]
    test.y=train.y1[(folds[[k]])]
    test.x=train.x1[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = gamma_list[i], probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])$auc
    mean.auc[i]=mean.auc[i]+r5/10
  }
}

m1=svm(x=as.matrix(train.x1),y=as.factor(train.y1), gamma = 0.1, probability = T)
class=predict(m1,test.x1,probability = T)
x=attributes(class)
r1=roc(test.y1~x$probabilities[,1], ci=T)

folds=createFolds(c(1:374), k=10)
gamma_list=seq(0.1, 1, 0.1)

mean.auc=rep(0,10)
for (i in 1:10){
  for (k in 1:10){
    train.y=train.y2[-(folds[[k]])]
    train.x=train.x2[-(folds[[k]]),]
    test.y=train.y2[(folds[[k]])]
    test.x=train.x2[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = gamma_list[i], probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])$auc
    mean.auc[i]=mean.auc[i]+r5/10
  }
}

m2=svm(x=as.matrix(train.x2),y=as.factor(train.y2), gamma = 0.1, probability = T)
class=predict(m2,test.x2,probability = T)
x=attributes(class)
r2=roc(test.y2~x$probabilities[,1], ci=T)


folds=createFolds(c(1:361), k=10)
gamma_list=seq(0.1, 1, 0.1)

mean.auc=rep(0,10)
for (i in 1:10){
  for (k in 1:10){
    train.y=train.y3[-(folds[[k]])]
    train.x=train.x3[-(folds[[k]]),]
    test.y=train.y3[(folds[[k]])]
    test.x=train.x3[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = gamma_list[i], probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])$auc
    mean.auc[i]=mean.auc[i]+r5/10
  }
}

m3=svm(x=as.matrix(train.x3),y=as.factor(train.y3), gamma = 0.1, probability = T)
class=predict(m3,test.x3,probability = T)
x=attributes(class)
r3=roc(test.y3~x$probabilities[,1], ci=T)


folds=createFolds(c(1:426), k=10)
gamma_list=seq(0.1, 1, 0.1)

mean.auc=rep(0,10)
for (i in 1:10){
  for (k in 1:10){
    train.y=train.y4[-(folds[[k]])]
    train.x=train.x4[-(folds[[k]]),]
    test.y=train.y4[(folds[[k]])]
    test.x=train.x4[(folds[[k]]),]
    m1=svm(x=as.matrix(train.x),y=as.factor(train.y), gamma = gamma_list[i], probability = T)
    class=predict(m1,test.x,probability = T)
    x=attributes(class)
    r5=roc(test.y~x$probabilities[,1])$auc
    mean.auc[i]=mean.auc[i]+r5/10
  }
}

m4=svm(x=as.matrix(train.x4),y=as.factor(train.y4), gamma = 0.1, probability = T)
class=predict(m4,test.x4,probability = T)
x=attributes(class)
r4=roc(test.y4~x$probabilities[,1], ci=T)


mean.auc=(r1$auc*109+r2$auc*152+r3$auc*165+r4$auc*100)/526
upper.auc=(r1$ci[3]*109+r3$ci[3]*152+r3$ci[3]*165+r4$ci[3]*100)/526
lower.auc=(r1$ci[1]*109+r3$ci[1]*152+r3$ci[1]*165+r4$ci[1]*100)/526

Check_density=function(v,x,y){
  x1=x[1:(length(x)-1)]
  x2=x[2:length(x)]
  return(max(y[(x1<=v)&(x2>=v)]))
}


c1=ci(r1, of="thresholds", thresholds="all")
c2=ci(r2, of="thresholds", thresholds="all")
c3=ci(r3, of="thresholds", thresholds="all")
c4=ci(r4, of="thresholds", thresholds="all")
upper1=(c1$sensitivity)[,3]
upper2=(c2$sensitivity)[,3]
upper3=(c3$sensitivity)[,3]
upper4=(c4$sensitivity)[,3]
lower1=(c1$sensitivity)[,1]
lower2=(c2$sensitivity)[,1]
lower3=(c3$sensitivity)[,1]
lower4=(c4$sensitivity)[,1]


specificity=unique(c(r1$specificities,r2$specificities,r3$specificities,r4$specificities))
specificity=sort(specificity)
sensitivity=upper=lower=rep(0,length(specificity))
for (i in 1:length(specificity)){
  m1=(109*Check_density(specificity[i],r1$specificities,r1$sensitivities)+
        152*Check_density(specificity[i],r2$specificities,r2$sensitivities)+
        165*Check_density(specificity[i],r3$specificities,r3$sensitivities)+
        100*Check_density(specificity[i],r4$specificities,r4$sensitivities))/526
  m2=(109*Check_density(specificity[i],r1$specificities,upper1)+
        152*Check_density(specificity[i],r2$specificities,upper2)+
        165*Check_density(specificity[i],r3$specificities,upper3)+
        100*Check_density(specificity[i],r4$specificities,upper4))/526
  m3=(109*Check_density(specificity[i],r1$specificities,lower1)+
        152*Check_density(specificity[i],r2$specificities,lower2)+
        165*Check_density(specificity[i],r3$specificities,lower3)+
        100*Check_density(specificity[i],r4$specificities,lower4))/526
  sensitivity[i]=m1
  upper[i]=m2
  lower[i]=m3
}

dat=data.frame(sensitivity,upper,lower,specificity=1-specificity)
dat=rbind(dat,c(0,0,0,0))

ggplot(aes(x=specificity,y=sensitivity),data=dat)+
  geom_path(size=1.1, color="#206999")+theme_bw()+theme(legend.text = element_text(size=14))+
  geom_ribbon(data=dat,aes(ymin=lower,ymax=upper),alpha=0.2)+
  geom_abline(slope=1, intercept=0, linetype=2)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=17,face="italic"))+
  theme(plot.title = element_text(hjust = 0.5),legend.position = c(0.7,0.15))+
  theme(legend.title = element_blank())+xlab("1-Specificity")+ylab("Sensitivity")









############################################################
##                                                        ##
##         Functional Analysis Part                       ##
##                                                        ##
############################################################


############################# Meta-analysis on GO genefamily ################################

feng=read.csv("feng.csv")               ### Read file in the "GO genefamily abundance"
vogtmann=read.csv("vogtmann.csv")
yu=read.csv("yu.csv")
zeller=read.csv("zeller.csv")


feng=feng[feng$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(feng$age,feng$BMI,feng$gender)
feng=feng[complete.cases(phenotype),]
feng$grouped_diagnosis=as.character(feng$grouped_diagnosis)
feng.GO=feng[,-c(1:14)]
x=(feng.GO>0)*1
sum.feng=colSums(x)



zeller=zeller[zeller$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(zeller$age,zeller$BMI,zeller$gender)
zeller=zeller[complete.cases(phenotype),]
zeller$grouped_diagnosis=as.character(zeller$grouped_diagnosis)
zeller.GO=zeller[,-c(1:14)]
x=(zeller.GO>0)*1
sum.zeller=colSums(x)



yu=yu[yu$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(yu$age,yu$BMI,yu$gender)
yu=yu[complete.cases(phenotype),]
yu$grouped_diagnosis=as.character(yu$grouped_diagnosis)
yu.GO=yu[,-c(1:14)]
x=(yu.GO>0)*1
sum.yu=colSums(x)



vogtmann=vogtmann[vogtmann$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(vogtmann$age,vogtmann$BMI,vogtmann$gender)
vogtmann=vogtmann[complete.cases(phenotype),]
vogtmann$grouped_diagnosis=as.character(vogtmann$grouped_diagnosis)
vogtmann.GO=vogtmann[,-c(1:14)]
x=(vogtmann.GO>0)*1
sum.vogtmann=colSums(x)


ratio=0.5
x1=(sum.feng>nrow(feng.GO)*ratio)*1
x2=(sum.vogtmann>nrow(vogtmann.GO)*ratio)*1
x3=(sum.yu>nrow(yu.GO)*ratio)*1
x4=(sum.zeller>nrow(zeller.GO)*ratio)*1
x=data.frame(x1,x2,x3,x4)
x_sum=4*apply(x, 1, mean)

feng.GO=feng.GO[,(x_sum==4)]
vogtmann.GO=vogtmann.GO[,(x_sum==4)]
yu.GO=yu.GO[,(x_sum==4)]
zeller.GO=zeller.GO[,(x_sum==4)]


## Austrian Cohort
response=as.factor(feng$grouped_diagnosis)
direction=p_value=p1=p2=rep(0,ncol(feng.GO))
for (i in 1:ncol(feng.GO)){
  m1=wilcox_test(feng.GO[,i]~response,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.feng=p1
p2.feng=p2


## American Cohort 
response=as.factor(vogtmann$grouped_diagnosis)
direction=p_value=p1=p2=rep(0,ncol(vogtmann.GO))
for (i in 1:ncol(vogtmann.GO)){
  m1=wilcox_test(vogtmann.GO[,i]~response,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.vogtmann=p1
p2.vogtmann=p2



## Chinese Cohort 
age_group=(yu$age %in% c(20:50))*1+(yu$age %in% c(51:60))*2+
  (yu$age %in% c(61:70))*3+(yu$age>70)*4
age_group=as.factor(age_group)
response=as.factor(yu$grouped_diagnosis)

direction=p_value=p1=p2=rep(0,ncol(yu.GO))
for (i in 1:ncol(yu.GO)){
  m1=wilcox_test(yu.GO[,i]~response|age_group,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.yu=p1
p2.yu=p2

## German Cohort 
age_group=(zeller$age %in% c(20:50))*1+(zeller$age %in% c(51:60))*2+
  (zeller$age %in% c(61:70))*3+(zeller$age>70)*4
age_group=as.factor(age_group)
response=as.factor(zeller$grouped_diagnosis)

direction=p_value=p1=p2=rep(0,ncol(zeller.GO))
for (i in 1:ncol(zeller.GO)){
  m1=wilcox_test(zeller.GO[,i]~response|age_group,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.zeller=p1
p2.zeller=p2


comb_p.increase=rep(0,ncol(feng.GO))
for (i in 1:ncol(feng.GO)){
  m1=wilkinsonp(c(p1.feng[i],p1.vogtmann[i],p1.zeller[i],p1.yu[i]),r=4)
  comb_p.increase[i]=m1$p
}
comb_p.increase=p.adjust(comb_p.increase,method = "BH")

comb_p.decrease=rep(0,ncol(feng.GO))
for (i in 1:ncol(feng.GO)){
  m1=wilkinsonp(c(p2.feng[i],p2.vogtmann[i],p2.zeller[i],p2.yu[i]),r=4)
  comb_p.decrease[i]=m1$p
}
comb_p.decrease=p.adjust(comb_p.decrease,method = "BH")

## generate the CRC enriched/depleted GO genefamilies
GO.enriched=names(feng.GO)[comb_p.increase<0.05]
GO.depleted=names(feng.GO)[comb_p.decrease<0.05]







############################# Meta-analysis on KO genefamily ################################

feng=read.csv("feng.csv")               ### Read file in the "KO genefamily abundance"
vogtmann=read.csv("vogtmann.csv")
yu=read.csv("yu.csv")
zeller=read.csv("zeller.csv")


feng=feng[feng$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(feng$age,feng$BMI,feng$gender)
feng=feng[complete.cases(phenotype),]
feng$grouped_diagnosis=as.character(feng$grouped_diagnosis)
feng.KO=feng[,-c(1:14)]
x=(feng.KO>0)*1
sum.feng=colSums(x)



zeller=zeller[zeller$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(zeller$age,zeller$BMI,zeller$gender)
zeller=zeller[complete.cases(phenotype),]
zeller$grouped_diagnosis=as.character(zeller$grouped_diagnosis)
zeller.KO=zeller[,-c(1:14)]
x=(zeller.KO>0)*1
sum.zeller=colSums(x)



yu=yu[yu$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(yu$age,yu$BMI,yu$gender)
yu=yu[complete.cases(phenotype),]
yu$grouped_diagnosis=as.character(yu$grouped_diagnosis)
yu.KO=yu[,-c(1:14)]
x=(yu.KO>0)*1
sum.yu=colSums(x)



vogtmann=vogtmann[vogtmann$grouped_diagnosis %in% c("CRC","normal_control"),]
phenotype=data.frame(vogtmann$age,vogtmann$BMI,vogtmann$gender)
vogtmann=vogtmann[complete.cases(phenotype),]
vogtmann$grouped_diagnosis=as.character(vogtmann$grouped_diagnosis)
vogtmann.KO=vogtmann[,-c(1:14)]
x=(vogtmann.KO>0)*1
sum.vogtmann=colSums(x)


ratio=0.5
x1=(sum.feng>nrow(feng.KO)*ratio)*1
x2=(sum.vogtmann>nrow(vogtmann.KO)*ratio)*1
x3=(sum.yu>nrow(yu.KO)*ratio)*1
x4=(sum.zeller>nrow(zeller.KO)*ratio)*1
x=data.frame(x1,x2,x3,x4)
x_sum=4*apply(x, 1, mean)

feng.KO=feng.KO[,(x_sum==4)]
vogtmann.KO=vogtmann.KO[,(x_sum==4)]
yu.KO=yu.KO[,(x_sum==4)]
zeller.KO=zeller.KO[,(x_sum==4)]


## Austrian Cohort
response=as.factor(feng$grouped_diagnosis)
direction=p_value=p1=p2=rep(0,ncol(feng.KO))
for (i in 1:ncol(feng.KO)){
  m1=wilcox_test(feng.KO[,i]~response,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.feng=p1
p2.feng=p2


## American Cohort 
response=as.factor(vogtmann$grouped_diagnosis)
direction=p_value=p1=p2=rep(0,ncol(vogtmann.KO))
for (i in 1:ncol(vogtmann.KO)){
  m1=wilcox_test(vogtmann.KO[,i]~response,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.vogtmann=p1
p2.vogtmann=p2



## Chinese Cohort 
age_group=(yu$age %in% c(20:50))*1+(yu$age %in% c(51:60))*2+
  (yu$age %in% c(61:70))*3+(yu$age>70)*4
age_group=as.factor(age_group)
response=as.factor(yu$grouped_diagnosis)

direction=p_value=p1=p2=rep(0,ncol(yu.KO))
for (i in 1:ncol(yu.KO)){
  m1=wilcox_test(yu.KO[,i]~response|age_group,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.yu=p1
p2.yu=p2

## German Cohort 
age_group=(zeller$age %in% c(20:50))*1+(zeller$age %in% c(51:60))*2+
  (zeller$age %in% c(61:70))*3+(zeller$age>70)*4
age_group=as.factor(age_group)
response=as.factor(zeller$grouped_diagnosis)

direction=p_value=p1=p2=rep(0,ncol(zeller.KO))
for (i in 1:ncol(zeller.KO)){
  m1=wilcox_test(zeller.KO[,i]~response|age_group,alternative="greater")
  if (pvalue(m1)<0.5) direction[i]=1
  else direction[i]=-1
  p1[i]=pvalue(m1)
  p2[i]=1-pvalue(m1)
  p_value[i]=2*min(pvalue(m1),1-pvalue(m1))
}
p1.zeller=p1
p2.zeller=p2


comb_p.increase=rep(0,ncol(feng.KO))
for (i in 1:ncol(feng.KO)){
  m1=wilkinsonp(c(p1.feng[i],p1.vogtmann[i],p1.zeller[i],p1.yu[i]),r=4)
  comb_p.increase[i]=m1$p
}
comb_p.increase=p.adjust(comb_p.increase,method = "BH")

comb_p.decrease=rep(0,ncol(feng.KO))
for (i in 1:ncol(feng.KO)){
  m1=wilkinsonp(c(p2.feng[i],p2.vogtmann[i],p2.zeller[i],p2.yu[i]),r=4)
  comb_p.decrease[i]=m1$p
}
comb_p.decrease=p.adjust(comb_p.decrease,method = "BH")

## generate the CRC enriched/depleted KO genefamilies
KO.enriched=names(feng.KO)[comb_p.increase<0.05]
KO.depleted=names(feng.KO)[comb_p.decrease<0.05]







