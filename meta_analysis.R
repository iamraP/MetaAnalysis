#install.packages((c"robumeta","metafor","dplyr")
#install.packages("rpsychi")

library(rpsychi)
library(robumeta)
library(metafor)
library(dyplr)
library(effectsize)
library(tidyr)


setwd("G:/Meine Ablage/PhD/FMT_MetaAnalysis")

data <- read.csv("test_export.csv",sep =",",fileEncoding="UTF-8-BOM")
data <- drop_na(data)
data$v <- as.numeric(data$v)
data_prep <- escalc(measure ="SMD", ,ni =N,yi=d, vi=v,data= data ,slab=study.ID)

res <- rma(yi=d,vi=v, data=data)
predict(res)
forest(res,cex=.8,header=c("Author(s), Year, [subgroup]"))




rma(data$d.,)

#create a string, with tau², I² and the Q-statistic, which can be added to the plot
rma_string <- function(rma_res) {
  text <-paste("Tau²: ",toString(round(rma_res$tau2,digits=4)),", (SE: ",toString(round(rma_res$se.tau2,digits=2)),")\nI²: ",toString(round(rma_res$I2,digits=2)),"%\nQ(df=",toString(round(rma_res$k-rma_res$p,digits=0)),")=",toString(round(rma_res$QE,digits=2)),", p=",toString(round(rma_res$QEp,digits=5)))
  return(text)
}

### Prepare data for analysis ##########
#1 way anova from mean and sample for deZambotti2012
deZamb <- read.csv("data/deZambotti2012.csv",sep =",",fileEncoding="UTF-8-BOM")
deZamb_anova <- with(deZamb, ind.oneway.second(Mean,SD,n))


## Load Data for Meta anlysis###

within <- read.csv("data/within_session.csv",sep =",",fileEncoding="UTF-8-BOM")
across <- read.csv("data/across_sessions.csv",sep =",",fileEncoding="UTF-8-BOM")


# calculate z correlation ( original variable should be expressed as the Fisher’s r-to-z transformed correlation coefficient)
within <- escalc(measure ="ZCOR",ri =r ,ni =n,data= within,slab=Name)
across <- escalc(measure ="ZCOR",ri =r ,ni =n,data= across,slab=Name)

within$multiple_bands_bin <- within$multiple_bands!="-"
across$multiple_bands_bin <- across$multiple_bands!="-"

### Within Session Modulation ################## 
# calculate model (random effects)
res_w <- rma(yi,vi, data=within)
res_w
predict(res_w) # confidence intervalls
confint(res_w)

#check if any of the studies is too influential?
inf <- influence(res_w)
print(inf)
plot(inf)

baujat(res_w)  

sav <-gosh(res_w, subset=20000)
plot(sav, out=3,breaks=100)
#create a string, with tau², I² and the Q-statistic, which can be added to the plot
report_w <- rma_string(res_w)

#convert z-to r add additonal info 
forest(res_w, atransf=transf.ztor,cex=.8,header=c("Author(s), Year, [subgroup]"),ilab=cbind(within$n, within$train_direction, within$multiple_bands,within$sessions), ilab.xpos=c(-9.4,-7.6,-5.5,-3.3),xlim=c(-16,7),alim=c(-3,2),ylim=c(-4,11))
text(c(-9.4,-7.6,-5.5,-3.3),c(10,10,10,10.05),cex=.8, c("n", "modulation\ndirection", "additonal\nbands", "sessions"))
text(c(-13.5),c(11),cex=1, c("Modulation within Sessions"))
text(c(-12),c(-2.7),cex=0.8,report_w )

#publication bias? 
funnel(res_w)
funnel(trimfill(res_w), las=1)

regtest(res_w)

#moderater analysis

res_w.moddir <- rma(yi,vi, mods= ~train_direction, data =within)
res_w.moddir

res_w.modses <- rma(yi,vi, mods = ~sessions, data = within)
res_w.modses

res_w.modban <- rma(yi,vi, mods= ~multiple_bands_bin, data =within)
res_w.modban

### Across Session Modulation ################## 
# calculate model (random effects)
res_a <- rma(yi,vi, data=across)
res_a
predict(res_a) # confidence intervalls
confint(res_a) # if confidence intervalls contain 0 values ther is no between study heterogeneity?

#check if any of the studies is too influential?
inf_a <- influence(res_a)
print(inf_a)
plot(inf_a)

#create a string, with tau², I² and the Q-statistic, which can be added to the plot
report_a <- rma_string(res_a)

#convert z-to r add additonal info 
forest(res_a, atransf=transf.ztor,cex=.8,header=c("Author(s), Year, [subgroup]"),ilab=cbind(across$n, across$train_direction, across$multiple_bands,across$sessions), ilab.xpos=c(-9.4,-7.6,-5.5,-3.3),xlim=c(-16,7),alim=c(-3,2),ylim=c(-4,10))
text(c(-9.4,-7.6,-5.5,-3.3),c(8.8,8.8,8.8,8.85),cex=.8, c("n", "modulation\ndirection", "additonal\nbands", "sessions"))
text(c(-13.1),c(10),cex=1, c("Modulation across Sessions"))
text(c(-12),c(-2.7),cex=0.8,report_a)


reporter(res_a, format="word")
#publication bias? 
funnel(res_a)
regtest(res_a)

#moderater analysis

res_a.moddir <- rma(yi,vi, mods= ~train_direction, data =across)
res_a.moddir

res_a.modses <- rma(yi,vi, mods = ~sessions, data = across)
res_a.modses

res_a.modban <- rma(yi,vi, mods= ~multiple_bands_bin, data =across)
res_a.modban

