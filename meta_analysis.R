install.packages(c("robumeta","metafor","dplyr"))
#install.packages("rpsychi")

library(rpsychi)
library(robumeta)
library(metafor)
library(dplyr)
library(effectsize)
library(tidyr)



### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, x) {
  bquote(paste(.(text),
                    " (Q = ", .(formatC(x$QE, digits=2, format="f")),
                    ", df = ", .(x$k - x$p),
                    ", p ", .(metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep=" ")), "; ",
                    "I²", " = ", .(formatC(x$I2, digits=1, format="f")), "%, ",
                    "tau²", " = ", .(formatC(x$tau2, digits=2, format="f")), ")"))}






#create a string, with tau², I² and the Q-statistic, which can be added to the plot
rma_string <- function(rma_res) {
  text <-paste("Tau²: ",toString(round(rma_res$tau2,digits=4)),", (SE: ",toString(round(rma_res$se.tau2,digits=2)),")\nI²: ",toString(round(rma_res$I2,digits=2)),"%\nQ(df=",toString(round(rma_res$k-rma_res$p,digits=0)),")=",toString(round(rma_res$QE,digits=2)),", p=",toString(round(rma_res$QEp,digits=5)))
  return(text)rma
}




setwd("G:/Meine Ablage/PhD/FMT_MetaAnalysis")

data <- read.csv("export_effect.csv",sep =",",fileEncoding="UTF-8-BOM")
data$ID <- seq.int(nrow(data))

# create the subesets of data acquired in a special way to rename it to merge it later
mean_SD <- subset(data,!is.na(M_SD))
group <- subset(data,!is.na(d_group))
no_correction <- subset(data,!is.na(d_cor_0))
within_subject_chen <- subset(data,!is.na(d_cor_chen))
within_subject_conservative <- subset(data,!is.na(d_cor_097))
within_subject_liberal  <- subset(data,!is.na(d_cor_074))

#rename the columns so they match and drop the empty columns

mean_SD <- mean_SD %>% rename("d" = "M_SD", "v" = "var_M_SD")
group <- group %>% rename("d" = "d_group", "v" = "var_group")
no_correction <- no_correction %>% rename("d" = "d_cor_0", "v" = "v_cor_0")
within_subject_chen <- within_subject_chen %>% rename("d" = "d_cor_chen", "v" = "v_cor_chen")
within_subject_conservative <- within_subject_conservative %>% rename("d" = "d_cor_097", "v" = "v_cor_097")
within_subject_liberal  <-  within_subject_liberal %>% rename("d" = "d_cor_074", "v" = "v_cor_074")


mean_SD <- mean_SD %>% select(-contains(c("_cor","_group")))
group <- group %>% select(-contains(c("_cor","_group","M_SD")))
within_subject_no_correction <- no_correction %>% select(-contains(c("_cor","_group","M_SD")))
within_subject_chen <- within_subject_chen %>% select(-contains(c("_cor","_group","M_SD")))
within_subject_conservative <- within_subject_conservative %>% select(-contains(c("_cor","_group","M_SD")))
within_subject_liberal <- within_subject_liberal %>% select(-contains(c("_cor","_group","M_SD")))


#merge the dataframes - create four different ones, for the different correction modes

stable <- rbind(mean_SD,group) 
no_correction <- rbind(stable,within_subject_no_correction)
chen <- rbind(stable,within_subject_chen)
conservative <- rbind(stable,within_subject_conservative)
liberal <- rbind(stable,within_subject_liberal)

### loop over the 4 dataframes

df_list = list(no_correction, chen, conservative, liberal)
df_names = list("no_correction", "chen", "conservative", "liberal")
for (i in 1:4 ){
  data_set = as.data.frame(df_list[i])
  data_set$names <- paste(data_set$Author,as.character(data_set$Year),sep=", ")
  V<- vcalc(vi=as.numeric(v), cluster =study.ID , time1 = t1, time2 = t2, grp1 = grp1 , grp2 = grp2, w1 = grp2_size, w2=grp2_size, data=data_set, phi = 0.9)
  res <- rma.mv(yi=d,V=V, random = ~ 1 | study.ID/effect_type, data=data_set,slab =names)
  sav <- robust(res, cluster = study.ID, adjust=TRUE,clubSandwich =TRUE)
  predict(sav)
  forest(sav, cex=.8,header=c("Author(s), Year"),rows = rev(seq(1,34,2)),
       order =data_set$names,
       ilab=cbind(data_set$Direction, data_set$Session, data_set$Subjects, data_set$effect_type), 
       ilab.xpos=c(-13,-11,-9,-7),
       ilab.pos = 4,
       xlim=c(-18,7),
       alim=c(-3,3),
       ylim=c(-4,36), 
       mlab = eval(mlabfun("RE Model", sav)))
  text(c(-12,-10,-8,-6),c(36,36,36,36.05),cex=.8, c("modulation\ndirection", "Sessions", "Subjects", ">Effect"))
  text(c(-16.5),c(36.6),cex=1, c(as.character(df_names[i])))
  }


influence(res)
res


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
inf <- influence(res)
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

