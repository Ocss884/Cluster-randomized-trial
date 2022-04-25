library(tidyverse)
library(gee)
library(car)
library(geepack)
library(ggpubr)

set.seed(283)

#Monte Carlo replication times 
replicate = 1000              

#Generating clusters
N0 = 3000
C = 160
nsize = round(runif(C,N0 / C * 0.6, N0 / C * 1.4))
N = sum(nsize)

#Propensity scores for clusters
ec = 0.3                                             

#Define cluster index and potential outcome 
cluster = numeric(N)      
c = 1
for (i in 1:N)
{
  if (i > sum(nsize[1:c])) c = c + 1 
  cluster[i] = c
}

#x : individual level covariate
#Y1 : treatment potential outcome
#Y0 : control potential outcome
x = (cluster / C) + runif(N, -1, 1)    
# Center x at its mean as required in the paper
x = x - mean(x)                                                                  
Y1 = x^3 + 2 * nsize[cluster] / (N / C)+rnorm(N, 0, 1)                                   
Y0 = (cluster / C) + x^2 + rnorm(N, 0, 1)                                              

#Calculate sum to prepare for cluster-total regressions
Ysum1 = numeric(C)
Ysum0 = numeric(C)
xsum = numeric(C)

for (c in 1:C)
{
  Ysum1[c] = sum(Y1[cluster == c])
  Ysum0[c] = sum(Y0[cluster == c])
  xsum[c] = sum(x[cluster == c])
}
#Calculate centered cluster size to prepare for cluster total regressions
nc = nsize - N / C

#Calculate causal effect
causal = mean(Y1-Y0)

tauIRp=0;tauTR2p=0;tauARp=0;tauOLSp=0
tauIRpols=0;tauTR2pols=0;tauARpols=0;tauOLSpols=0
tauIRphw=0;tauOLSphw=0

#Store the estimators for each Monte Carlo

tauIR=numeric(replicate)
tauTR2=numeric(replicate)
tauAR=numeric(replicate)
tauOLS=numeric(replicate)

#Store the variance estimators for each Monte Carlo
tauIRev=numeric(replicate);   tauIRevols=numeric(replicate);  tauIRevhw=numeric(replicate)
tauTR2ev=numeric(replicate);  tauTR2evols=numeric(replicate)
tauARev=numeric(replicate);   tauARevols=numeric(replicate)
tauOLSev=numeric(replicate);  tauOLSevols=numeric(replicate); tauOLSevhw=numeric(replicate)

Z=numeric(N)
for (i in 1:replicate)
{
  #Generate treatment assignment and observed data
  Zc=numeric(C)
  Z0=sample(C,C*ec,replace=FALSE)
  for (j in 1:(C*ec)) Zc[Z0[j]]=1
  for (j in 1:N) Z[j]=Zc[cluster[j]]
  Y=Z*Y1+(1-Z)*Y0
  Ysum=Zc*Ysum1+(1-Zc)*Ysum0
  #ancova
  dat=data.frame(Y,Z,x,cluster)
  r=geeglm(Y~1+Z+x,id=cluster,corstr="independence",data=dat)
  tauOLS[i]=r$coefficients[2]
  tauOLSev[i]=summary(r)$coefficients[2,2]^2
  r=lm(Y~1+Z+x)
  tauOLSevols[i]=summary(r)$coefficients[2,2]^2  
  tauOLSevhw[i]=hccm(r,type="hc0")[2,2]    
  
  if ((tauOLS[i]-sqrt(tauOLSev[i])*qnorm(0.975)<causal)&(causal<tauOLS[i]+sqrt(tauOLSev[i])*qnorm(0.975))) tauOLSp=tauOLSp+1
  if ((tauOLS[i]-sqrt(tauOLSevols[i])*qnorm(0.975)<causal)&(causal<tauOLS[i]+sqrt(tauOLSevols[i])*qnorm(0.975))) tauOLSpols=tauOLSpols+1
  if ((tauOLS[i]-sqrt(tauOLSevhw[i])*qnorm(0.975)<causal)&(causal<tauOLS[i]+sqrt(tauOLSevhw[i])*qnorm(0.975))) tauOLSphw=tauOLSphw+1
  
  
  #Individual level regression with covariate x_{ij}
  dat=data.frame(Y,Z,x,cluster)
  r=geeglm(Y~1+Z+x+Z*x,id=cluster,corstr="independence",data=dat)
  tauIR[i]=r$coefficients[2]
  tauIRev[i]=(summary(r)$coefficients[2,2])^2
  r=lm(Y~1+Z+x+Z*x)
  tauIRevols[i]=summary(r)$coefficients[2,2]^2  
  tauIRevhw[i]=hccm(r,type="hc0")[2,2]  
  if ((tauIR[i]-sqrt(tauIRev[i])*qnorm(0.975)<causal)&(causal<tauIR[i]+sqrt(tauIRev[i])*qnorm(0.975))) tauIRp=tauIRp+1
  if ((tauIR[i]-sqrt(tauIRevols[i])*qnorm(0.975)<causal)&(causal<tauIR[i]+sqrt(tauIRevols[i])*qnorm(0.975))) tauIRpols=tauIRpols+1
  if ((tauIR[i]-sqrt(tauIRevhw[i])*qnorm(0.975)<causal)&(causal<tauIR[i]+sqrt(tauIRevhw[i])*qnorm(0.975))) tauIRphw=tauIRphw+1
  
  xm=xsum[cluster]/nsize[cluster]
  
  #Cluster total regression with cluster size and additional covariate
  r=lm(Ysum/(N/C)~1+Zc+nc+xsum+Zc*nc+Zc*xsum)   
  tauTR2[i]=r$coefficients[2]
  tauTR2ev[i]=hccm(r,type="hc0")[2,2]
  tauTR2evols[i]=summary(r)$coefficients[2,2]^2  
  if ((tauTR2[i]-sqrt(tauTR2ev[i])*qnorm(0.975)<causal)&(causal<tauTR2[i]+sqrt(tauTR2ev[i])*qnorm(0.975))) tauTR2p=tauTR2p+1
  if ((tauTR2[i]-sqrt(tauTR2evols[i])*qnorm(0.975)<causal)&(causal<tauTR2[i]+sqrt(tauTR2evols[i])*qnorm(0.975))) tauTR2pols=tauTR2pols+1
  
  #Average with additional covariate
  xmeanc=xsum/nsize-mean(xsum/nsize)
  r=lm(Ysum/nsize~1+Zc+xmeanc+Zc*xmeanc)   
  tauAR[i]=r$coefficients[2]
  tauARev[i]=hccm(r,type="hc0")[2,2] 
  tauARevols[i]=summary(r)$coefficients[2,2]^2
  if ((tauAR[i]-sqrt(tauARev[i])*qnorm(0.975)<causal)&(causal<tauAR[i]+sqrt(tauARev[i])*qnorm(0.975))) tauARp=tauARp+1
  if ((tauAR[i]-sqrt(tauARevols[i])*qnorm(0.975)<causal)&(causal<tauAR[i]+sqrt(tauARevols[i])*qnorm(0.975))) tauARpols=tauARpols+1
  
}

#Arrange the data into dataframe

datIR=rbind(mean(tauIR)-causal,sd(tauIR),mean(sqrt(tauIRevols)),mean(sqrt(tauIRevhw)),mean(sqrt(tauIRev)),
            sqrt(mean((tauIR-causal)^2)),tauIRpols,tauIRphw,tauIRp)

datTR2=rbind(mean(tauTR2)-causal,sd(tauTR2),mean(sqrt(tauTR2evols)),mean(sqrt(tauTR2ev)),
             sqrt(mean((tauTR2-causal)^2)),tauTR2pols,tauTR2p)

datAR=rbind(mean(tauAR)-causal,sd(tauAR),mean(sqrt(tauARevols)),mean(sqrt(tauARev)),
            sqrt(mean((tauAR-causal)^2)),tauARpols,tauARp)

datOLS=rbind(mean(tauOLS)-causal,sd(tauOLS),mean(sqrt(tauOLSevols)),mean(sqrt(tauOLSevhw)),mean(sqrt(tauOLSev)),
             sqrt(mean((tauOLS-causal)^2)),tauOLSpols,tauOLSphw,tauOLSp)

mul=100
#First row of output
dat1=round(mul*cbind(datIR,datOLS),2)
#Second row of output
dat2=round(mul*cbind(datTR2,datAR),2)

#Scale the dataframe

dat1[7,]=dat1[7,]/(mul*replicate)
dat1[8,]=dat1[8,]/(mul*replicate)
dat1[9,]=dat1[9,]/(mul*replicate)
dat2[6,]=dat2[6,]/(mul*replicate)
dat2[7,]=dat2[7,]/(mul*replicate)

#Plot the coverage rates
bar_data1 = as.data.frame(t(dat1[7:9,1:2]))
names(bar_data1) <- c("OLS", "HW", "LZ")

bar_data2 = as.data.frame(t(dat2[6:7,1:2]))
bar_data2$LZ <- rep(0, 2)
names(bar_data2) <- c("OLS", "HW", "LZ")

bar_data <- rbind(bar_data1, bar_data2)

bar_data <- rownames_to_column(bar_data, var = "row_name")

bar_data <- bar_data %>% pivot_longer(OLS:LZ) 

bar_data$name <- factor(bar_data$name, 
                        levels = c("OLS","HW", "LZ"))

bar_data$row_name <- factor(bar_data$row_name, 
                            levels = c("1", "2", "3", "4"))

rate <- ggplot(data = bar_data, aes(x = row_name, y = value, fill = name)) + 
  geom_bar(position="dodge", stat = "identity") + ylim(0, 1) + 
  geom_abline(slope = 0, intercept = 0.95, lty = 2) +
  scale_x_discrete(labels = c(expression(hat(tau)[I]),
                              expression({ANOVA}),
                              expression({hat(tau)[T]}),
                              expression(hat(tau)[pi]))) +
  theme(axis.text.x = element_text(angle = 60, size = 13)) + 
  labs(x = "", y = "Coverage rate")

rate
