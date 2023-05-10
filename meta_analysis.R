#library(rpsychi)
library(robumeta)
library(metafor)
library(dplyr)
library(effectsize)
library(tidyr)
library(ggplot2)
library(plotrix)


#### Prerequisites #######
####### new graphic window 
#dev.new()
par(xpd=NA)


###  helper functions to add Q-test

#reporting of univariate model heteorgenity
uni_hetero <- function(text, x) {
  ital_p <- ","~italic("p ")
  bquote(paste(.(text), " (Q(", .(x$k - x$p),") = ", .(formatC(x$QE, digits=2, format="f")),
               ", p ", .(metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep="")),
               ", I²", " = ", .(formatC(x$I2, digits=1, format="f")), "%, ",
               "τ²", " = ", .(formatC(x$tau2, digits=2, format="f")), ")"))}

#reporting of mutlivariate model heteorgenity
mv_hetero <- function(text, x, effects) {
  #calculate the I-sqaured 
  
  W <- diag(1/res$vi)
  X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  x$I2 <- 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
  
  sigma_strings <- c()
  for (i in 1:length(effects)){
   sigma_strings<- c(sigma_strings, paste(", σ²_", effects[i]," = ", (formatC(x$sigma2[i], digits=2, format="f")),sep=""))
  }
  sigma_string <- paste(sigma_strings,collapse="")
  
  p_value <- metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep=" ")
  Q_value <- formatC(x$QE, digits=2, format="f")
  
  paste(text, " (Q(", (x$k - x$p),") = ", (formatC(x$QE, digits=2, format="f")),
               ", p ", (metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep="")),
               ", I²", " = ", (formatC(x$I2, digits=1, format="f")), "%",
               sigma_string, ")")}

#color forest plot

color_forest_plot <- function(information, rows, x_pos,row_idx) {

  ##color coding for the graphs: 
  transperent_color_set_DIAD <- c("#5EB5188C","#B2E8878C","#E07B7D8C", "#B5181B8C" )
  transperent_color_set_dir <- c("#0099998C","#93ffff8C")
  transperent_color_set_sub <- c("#fec44f8C","#a750178C")
  transperent_color_set_ses <-  c("#99d8c98C","#2ca25f8C")
  transperent_color_set_freq <- c("#e34a338C", "#fdbb848C")
  transperent_color_set_af <- c("#9ebcda8C", "#8856a78C")
  

  j=1
  #direction
  dir_col <- match(information[j],c("↑","↓"))
  rect(x_pos[j]-0.25, rows[row_idx]-0.46, x_pos[j]+0.25, rows[row_idx]+0.46,density = NA, col = transperent_color_set_dir[dir_col])
  #effect subjects
  j=2
  sub_col <- match(information[j],c("W","A"))
  if (!is.na(sub_col)){
    rect(x_pos[j]-0.25, rows[row_idx]-0.46, x_pos[j]+0.68, rows[row_idx]+0.46,density = NA, col = transperent_color_set_sub[sub_col])
  }
  #effect session
  j=4
  ses_col <- match(information[j],c("W","A"))
  if (!is.na(ses_col)){
    rect(x_pos[j]-0.25, rows[row_idx]-0.46, x_pos[j]+0.68, rows[row_idx]+0.46,density = NA, col = transperent_color_set_ses[ses_col])
  }
  #frequency
  j=6
  freq_col <- match(information[j],c("F","I"))
  if (!is.na(freq_col)){
    rect(x_pos[j]-0.22, rows[row_idx]-0.46, x_pos[j]+0.22, rows[row_idx]+0.46,density = NA, col = transperent_color_set_freq[freq_col])
  }
  #add frequency
  j=7
  af_col <- match(information[j],c("α","β"))
  if (!is.na(af_col)){
    rect(x_pos[j]-0.22, rows[row_idx]-0.46, x_pos[j]+0.22, rows[row_idx]+0.46,density = NA, col = transperent_color_set_af[af_col])
  }
  #DIAD
  for (j in 8:length(information) ){
    diad_col <- match(information[j],c("++","+","-","--"))
    if (!is.na(diad_col)){
      rect(x_pos[j]-0.25, rows[row_idx]-0.46, x_pos[j]+0.25, rows[row_idx]+0.46,density = NA, col = transperent_color_set_DIAD[diad_col])
    }
  }
}
  




############# SETTINGS #############
correction_mode = 3
remove_outliers <- TRUE

effects <- list("mod", "bsl_pre_post", "bsl_mod")
correction <- c("chen", "con", "lib")


##Get Data 


setwd("G:/Meine Ablage/PhD/Neurofeedback/FMT_MetaAnalysis")
data <- read.csv("export_effect.csv",sep =",",fileEncoding="UTF-8-BOM")


#### Start Analysis

# index data, so that effects can be disentangled later on
data$index <- 1:nrow(data)


# #relable data 
data[data=="Across"] <- "A"
data[data=="Within"] <- "W"   
data$Direction[data$Direction==-1] <- "↓"
data$Direction[data$Direction==1] <- "↑"
data$IF[data$IF==0] <- "F"
data$IF[data$IF==1] <- "I"
data$add_freq[data$add_freq=="alpha"] <- "α"
data$add_freq[data$add_freq=="beta"] <- "β"
data[data=="Y"] <- "++"
data[data=="mY"] <- "+"
data[data=="mN"] <- "-"
data[data=="N"] <- "--"




print(paste("############## ", correction[correction_mode]," ###################### "))

# select values corresponding to the correction and rename them uniformly
data_set <- data %>% select(-contains(correction[-correction_mode]))
hedges_g <- paste("g",correction[correction_mode],sep="_")
variance <-paste("var_g",correction[correction_mode],sep="_")
data_set <- data_set %>% rename("g" = hedges_g, "v" = variance)

# only keep important data
data_set <- data_set[c("Author",
                       "index",
                       "Year",
                       "grp1","grp1_size",
                       "grp2","grp2_size",
                       "t1","t2",
                       "study.ID",
                       "n_freq",
                       "IF",
                       "g",
                       "v",
                       "effect_type",
                       "Direction",
                       "Subjects",
                       "add_freq",
                       "Session",
                       "number.of.sessions",
                       "PaperID",
                       "study.ID",
                       "A",
                       "B",
                       "C",
                       "D")]

# add name tags for the plot
data_set$names <- paste(data_set$Author,as.character(data_set$Year),sep=", ")

# drop studies without data, i.e. those w/o effect sizes -> douible in the figure if the ones listed at "not included in the model" are just the ones without effects sizes
data_set <- drop_na(data_set)

## Overall Model ##########

data_set <- data_set[order(data_set$effect_type),] #data_set[order(as.character(data_set$effect_type))]
#Variance
V<- vcalc(vi=as.numeric(v),
          cluster =study.ID ,
          time1 = t1,
          time2 = t2,
          grp1 = grp1,
          grp2 = grp2,
          w1 = grp2_size,
          w2 = grp2_size,
          data = data_set,
          phi = 0.9)

#the model
res <- rma.mv(yi = as.numeric(g),
              V = V,
              random = ~ 1 | study.ID/effect_type,
              data=data_set,
              slab =names)

assign(paste("res_rob",correction[correction_mode],sep = "."), res)

#the cluster robust estimation
res_rob <- robust(res,
                  cluster = study.ID,
                  adjust = TRUE,
                  clubSandwich = TRUE)

assign(paste("res",correction[correction_mode],sep = "."), res)

## Subgroup Models  #########

subsets <- list(### Subeffects
  subset(data_set,effect_type =="mod"),
  subset(data_set,effect_type =="bsl_pre_post"),
  subset(data_set,effect_type =="bsl_mod"),    #incl. outlier
  #no of sessions
  subset(data_set,number.of.sessions == 1),  #single session
  subset(data_set,number.of.sessions > 1),  #multiple   #incl. outlier
  # Modulation direction
  subset(data_set,Direction == "↓"),         # Down
  subset(data_set,Direction == "↑"),             #Up   #incl. outlier
  # individual frequencies             #
  subset(data_set,IF == "I"),
  subset(data_set,IF == "F"),             # standard freq    #incl. outlier
  # number of frequencies,    
  subset(data_set,n_freq > 1),         # Muti Band
  subset(data_set,n_freq == 1))                  # Single Band #incl. outlier      


# The order needs to match the order of subsets!!  maybe there's a way to implement this as a dict
subsets_label <- c("mod",
                   "bsl_pre_post",
                   "bsl_mod",
                   "single_sess",
                   "multi_sess",
                   "down",
                   "up",
                   "ind_freq",
                   "fix_freq",
                   "multi_band",
                   "single_band")




# counter for the amount of rows needed in the forest plot
n_grp_rows <- c()

outlier_subsets <- list() # subsets without influential outliers - Inspect all of them manually before excluding them!
outlier_subsets_labels <-list()  # labels for subsets without influential outliers


# iterate over all subsets and control for influential outliers#
#save the results in a text file
sink(file = "inf_out_dia.txt")
for (k in 1:length(subsets)){
  # print name of subset
  sub_set <-subsets[[k]]
  cat("______________________________________________________________________________________\n______________________________________________________________________________________\n                              ",subsets_label[k], "\n______________________________________________________________________________________\n")
  print(subsets_label[k])
  
  # remember the number of rows for the first three subsets, as they are inclued directly in the forest plot
  if (k<=3){
    n_grp_rows <- c(n_grp_rows,nrow(sub_set))
  }
  
  #calculate a multivariate model whenever more than one effect is derived from a single study or a the same type of subeffect is included, if the subgroup is not defined by the subeffect
  if ((sum(duplicated(sub_set$study.ID)) > 0 | (sum(duplicated(sub_set$effect_type)) >0 & !subsets_label[k] %in% effects))) {
    V<- vcalc(vi=as.numeric(v),
              cluster =study.ID ,
              time1 = t1,
              time2 = t2,
              grp1 = grp1,
              grp2 = grp2,
              w1 = grp2_size,
              w2 = grp2_size,
              data = sub_set,
              phi = 0.9)
    #set random effects according to included  studies
    if ((sum(duplicated(sub_set$effect_type)) >0 & !subsets_label[k] %in% effects) & sum(duplicated(sub_set$study.ID)) > 0 ){ 
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | study.ID/effect_type,
                    data=sub_set,
                    slab =names)
    }else if (sum(duplicated(sub_set$study.ID)) > 0 ){
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | study.ID,
                    data=sub_set,
                    slab =names)
    }else if (sum(duplicated(sub_set$effect_type)) >0 & !subsets_label[k] %in% effects){
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | effect_type,
                    data=sub_set,
                    slab =names)
    }
    # create cluster robust estimator
    res <- robust(res,
                  cluster = study.ID,
                  adjust = TRUE,
                  clubSandwich = TRUE)
    
    # calculate influence and outlier diganostics
    cook <-cooks.distance.rma.mv(res)
    hat <- hatvalues.rma.mv(res)
    res_non_robust <-rma(yi = as.numeric(g),
                         vi=as.numeric(v),
                         data=sub_set,
                         slab =names)
    
  # use a univariate model if none of the above applies 
  } else {
    res <- rma(yi = as.numeric(g),
               vi=as.numeric(v),
               data=sub_set,
               slab =names)
    cook <-cooks.distance.rma.uni(res)
    hat <- hatvalues.rma.uni(res)
    res_non_robust <- res
  }
  

  #calulate the Students residuals  (with non-robust model)
  stud_res <- residuals.rma(res_non_robust,type="rstudent")
  
  #save all models named res.correction_mode.subset
  assign(paste("res",correction[correction_mode],subsets_label[k],sep = "."), res)
  
  
  #save outliers and influencers 
  outliers <- names(stud_res[sqrt(stud_res^2)>1.96])
  influencers <- names(cook[sqrt(cook^2)>median(cook)+6*IQR(cook)])         
  
  #print the outliers and influencers to the text file
  print(paste("cutoff for outliers: 1.96"))
  print(paste("outliers: ", outliers))
  

  print(paste("cutoff for influencers: ", median(cook)+6*IQR(cook)))
  print(paste("influencers: ", influencers))
  
  # if outliers and influencer in both groups add them to list  to recalculate the model without them
  if (any(outliers %in% influencers)){
    outlier_subsets <-  append(outlier_subsets,list(subset(sub_set,!(sub_set$names %in% outliers[outliers %in% influencers]))))
    outlier_subsets_labels <-  append(outlier_subsets_labels, subsets_label[k])
  }
}

sink(file = NULL)


# ############### SUBSET OUTLIER EVALUATION - INVESTIGATE BY HAND!! ############## 

# res_investigate <- res_rob
# 
# # # estimate variance
# profile(res_investigate,sigma=1)
# dev.print(pdf, file=eval(paste("Diagnostics/profile_studyID_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))
# profile(res_investigate,sigma=2)
# dev.print(pdf, file=eval(paste("Diagnostics/profile_effectType_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))
# 
# 
# 
# 
# 

######## same analysis w/o the influential outlier #####

# Overall Model 

data_set_cleaned <- subset(data_set,data_set$study.ID!="RJHQE9FU") #make this more flexible
print("Influential Outliers removed!")

# calculate the overall effect 
#Variance
V_cleaned<- vcalc(vi=as.numeric(v),
                  cluster =study.ID ,
                  time1 = t1,
                  time2 = t2,
                  grp1 = grp1,
                  grp2 = grp2,
                  w1 = grp2_size,
                  w2 = grp2_size,
                  data = data_set_cleaned,
                  phi = 0.9)
#the model
res_cleaned <- rma.mv(yi = as.numeric(g),
                      V = V_cleaned,
                      random = ~ 1 | study.ID/effect_type,
                      data=data_set_cleaned,
                      slab =names)

#the cluster robust estimation
res_rob_cleaned <- robust(res_cleaned,
                          cluster = study.ID,
                          adjust = TRUE,
                          clubSandwich = TRUE)


# Subsets

for (k in 1:length(outlier_subsets)){
  sub_set <- outlier_subsets[[k]]
  if ((sum(duplicated(sub_set$study.ID)) > 0 | (sum(duplicated(sub_set$effect_type)) >0 & !outlier_subsets_labels[k] %in% effects) )) {
    V<- vcalc(vi=as.numeric(v),
              cluster =study.ID ,
              time1 = t1,
              time2 = t2,
              grp1 = grp1,
              grp2 = grp2,
              w1 = grp2_size,
              w2 = grp2_size,
              data = sub_set,
              phi = 0.9)
    #set random effects according to included  studies
    if (sum(duplicated(sub_set$effect_type) >0) & sum(duplicated(sub_set$study.ID)) > 0 & !outlier_subsets_labels[k] %in% effects ){ 
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | study.ID/effect_type,
                    data=sub_set,
                    slab =names)
    }else if (sum(duplicated(sub_set$study.ID)) > 0 ){
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | study.ID,
                    data=sub_set,
                    slab =names)
    }else if (sum(duplicated(sub_set$study.ID)) > 0 & !outlier_subsets_labels[k] %in% effects){
      res <- rma.mv(yi = as.numeric(g),
                    V = V,
                    random = ~ 1 | effect_type,
                    data=sub_set,
                    slab =names)
    }
    # create cluster robust estimator
    res <- robust(res,
                  cluster = study.ID,
                  adjust = TRUE,
                  clubSandwich = TRUE)
    # use a univariate model if none of the above applies 
  } else {
    res <- rma(yi = as.numeric(g),
               vi=as.numeric(v),
               data=sub_set,
               slab =names)
  }
  
  assign(paste("res_cleaned",correction[correction_mode],outlier_subsets_labels[k],sep = "."), res)
  
}



### Preparation for the  forest plot ####

# manually adjusted x positions for the labels, and the information 
x_pos <- c(-7.7,-6.9,-6.4,-5.6,-5.1,-4.3,-3.8,seq(7.74,9.9,0.63))


# calculate rows required in the forest plot for studies not included in the model 
n_grp_rows <- c(n_grp_rows,nrow(data) - sum(n_grp_rows))
# how many rows in total 
req_rows <- length(n_grp_rows)*2 + (length(n_grp_rows)-1)*2+ sum(n_grp_rows)

y_lim=c(-24, req_rows)

# row distribution
poly_rows<- c()
header_rows <- c(req_rows-3)
add_rows <- c()
forest_rows <- c()
stop <- req_rows-4
for (t in 1:length(n_grp_rows)){
  start <-stop-n_grp_rows[t]+1
  rows <- start:stop
  if (t < length(n_grp_rows)){
    forest_rows <- c(rows,forest_rows)
    stop <- start - 4
    header_rows <- c(stop+1, header_rows)
    poly_rows <- c(stop+2.8, poly_rows)
  } else { 
    header_rows[1] <- header_rows[1]-1
    add_rows<- c(rows-1,add_rows)
    
  }
}
# dev.new(width=1130, height=680,noRStudioGD = TRUE, unit="px")

# set the current model including all studies as the global model
global_model <- eval(parse(text=paste("res_rob",correction[correction_mode],sep = ".")))


### Plotting ####

#mar:  bottom, left, top and right
par(oma =c(0,0,0,0), mar= c(6,0,0,12),font=1,cex=0.8)

#create the actual plot 
waldi <-forest(global_model,
               xlim = c(-16.2,6.9),
               addpred = TRUE,
               col="#9e0000",
               ilab = cbind(data_set$Direction,
                            data_set$Subjects,
                            data_set$grp1_size,
                            data_set$Session,
                            data_set$number.of.sessions,
                            data_set$IF,
                            data_set$add_freq,
                            data_set$A,
                            data_set$B,
                            data_set$C,
                            data_set$D), 
               xlab = c(""),
               showweights = TRUE,
               ilab.xpos = x_pos,
               ylim=y_lim,
               steps = 9,
               at = seq(-4,4),
               order=effect_type,
               rows=forest_rows,
               mlab= eval(mv_hetero("REM for global effect",global_model,c("s","e"))),
               psize=0.9,
               header="Author(s) and Year")


#forest_order <- (order(global_model$data$effect_type))
for (idx in 1:length(forest_rows)){ 
  study <- global_model$data[idx,]
  information <- c(study$Direction,
                   study$Subjects,
                   study$grp1_size,
                   study$Session,
                   study$number.of.sessions,
                   study$IF,
                   study$add_freq,
                   study$A,
                   study$B,
                   study$C,
                   study$D)
  color_forest_plot(information, forest_rows, x_pos,idx)
}
  


#add missing studies
additional_studies <- data$index[!data$index %in% data_set$index]
for (study in 1:length(additional_studies)){ 
  add_study <- subset(data, data$index == additional_studies[study])
  information <- c(add_study$Direction,
                   add_study$Subjects,
                   add_study$grp1_size,
                   add_study$Session,
                   add_study$number.of.sessions,
                   add_study$IF,
                   add_study$add_freq,
                   add_study$A,
                   add_study$B,
                   add_study$C,
                   add_study$D)
  text(-9*1.8,add_rows[study],pos=4, paste(add_study$Author,add_study$Year,sep=", "),cex=0.8)
  for (j in 1:length(information) ){
    text(x_pos[j],add_rows[study],information[j],cex=0.8)
         }
  color_forest_plot(information, add_rows, x_pos, study)
}

  

#w/o outlier
addpoly(res_rob_cleaned, row=-2.4, mlab=eval(mv_hetero("REM w/o influential outliers",res_rob_cleaned,c("s","e"))),addpred=TRUE,col="#d28888")

# subgroup plot below 
y_pos_subgroups <- -2
y_pos_clean_subgroups <-c()
model_options <- c("res","res_cleaned")


for (i in 1: length(subsets_label)){
  # make sure the first three models are placed correctly
  if (i<=3){y_pos_subgroups <- poly_rows[4-i]} else if(i==4){y_pos_subgroups<- -2.1}
  
  # print outlier model immediatly after the non cleaned model
  if (subsets_label[i] %in% outlier_subsets_labels){ 
    models_to_print <- 2
  } else models_to_print <- 1 
  # leave space between different influence factors 
  if (i%%2==0 & i >3 ){y_pos_subgroups <- y_pos_subgroups -1} 
  for (k in 1:models_to_print){
    #write in the next row 
    if (i>3 | k>1) {y_pos_subgroups <- y_pos_subgroups -1.5}
    #get next subgroup
    subgroup = eval(parse(text=paste(model_options[k],correction[correction_mode],subsets_label[i],sep = ".")))
    #color cleaned models accordingly 
    if(k==2){col<-"#c4d2ee"
      description <- "REM w/o influential outliers for"
    } else {col <-"#013394"
      description <- "REM for"}
    if (is.null(subgroup$random)){
      addpoly(subgroup, y_pos_subgroups, mlab=eval(uni_hetero(paste(description, subsets_label[i], sep=" "),subgroup)),addpred=TRUE,col=col)
    } else if (grepl("effect_type",subgroup$random)) {
      addpoly(subgroup, y_pos_subgroups, mlab=eval(mv_hetero(paste(description, subsets_label[i], sep=" "),subgroup,c("s", "e"))),addpred=TRUE,col=col)
    } else{
      addpoly(subgroup, y_pos_subgroups, mlab=eval(mv_hetero(paste(description, subsets_label[i], sep=" "),subgroup,c("s"))),addpred=TRUE,col=col)
    }
  }
}


# Add headers for columns
text(c(-7.8,-6.6,-5.3,-4,8.5),33, pos=3,cex=.8, c("Direction",
                                                  "Subjects",
                                                  "Sessions",
                                                  "Frequency",
                                                  "Study DIAD - Risk of Bias"))
text(4,32, pos=3,cex=.8, c("Weight"))

# Add sub headers for columns
text(c(x_pos[2:6],-3.65 ,x_pos[8:12]),32,pos=3,cex=.7, c("Effect",
                                                         "n",
                                                         "Effect",
                                                         "n",
                                                         "F/I",
                                                         "Additional",
                                                         "A","B","C","D"))

### switch to bold italic font
par(font=2,cex=0.7)

### add text for the subgroups
text(-16, rev(header_rows), pos=4, c("Beginning vs. Ending of Neurofeedback",
                                     "Resting State Baseline vs. After Neurofeedback",
                                     "Resting State vs. Neurofeedback",
                                     "Not included in the model"))

par(xpd=NA)

# xlabel 
text(0,-30, "Observed Outcome (Hedge's g)",font=1)

#seperate study DIAD 
segments(6.85, -30, x1 = 6.85, y1 = 35,lty = "dotdash")

# add description for Study DIAD
text(c(rep(7.8,5)),seq(-2,-6),pos =4, c(
  "Ratings:",
  "++ : \"Yes\"",
  "+  : \"Maybe Yes\"",
  "-  : \"Maybe No\"",
  "-- : \"No\""),
  cex=0.9,font=1)




# 
# mark influential studies manually (currently set to scaling of A4)##### 
if (correction_mode==1){
  ### #  Chen #
  #text(-14.1,8.2,cex=0.8, c("†")) #Chen
  text(-5.6,5.5,cex=0.8, c("†"))  # bsl mod
  text(-5.2,-6,cex=0.8, c("†"))#single
  #Wang
  #text(-13.95,16.4,cex=1.4, c("*")) #Wang
  text(-5.34,-2,cex=1.4, c("*")) # overall model
  text(-5.42,5.7,cex=1.4, c("*"))  # bsl mod
  text(-5.35,-8.6,cex=1.4, c("*")) #multi
  text(-4.7,-14.2,cex=1.4, c("*")) #up
  text(-4.14,-19.7,cex=1.4, c("*")) #fix
  
  } else if (correction_mode == 2){

  ### #  Con #
  #text(-14.1,8.2,cex=0.8, c("†")) #Tseng
  text(-5.3,-20.7,cex=0.8, c("†"))#  multi
  #Wang
  #text(-13.95,16.4,cex=1.4, c("*")) #Wang
  text(-5.3,-2,cex=1.4, c("*")) #overall
  text(-5.4,-16.6,cex=1.4, c("*")) #fix
  text(-6.1,-11.2,cex=1.4, c("*")) #up
  
  }else if (correction_mode == 3){
  
  ### #  Lib  ##
  #text(-14.1,8.2,cex=0.8, c("†")) #Chen
  text(-5.6,5.5,cex=0.8, c("†"))  # bsl mod
  text(-5.25,-6,cex=0.8, c("†")) #single sess
  #Wang
  #text(-13.95,16.4,cex=1.4, c("*")) #Wang
  text(-5.34,-2,cex=1.4, c("*")) # overall model
  text(-5.4,5.7,cex=1.4, c("*"))  # bsl mod
  text(-5.45,-16.7,cex=1.4, c("*")) #fix
  text(-3.8,-22.2,cex=1.4, c("*")) #single band
}


# ### Assessment of Publication bias ######
# #Publication bias Funnel plot

#labels to numbers 
slab_funnel = data.frame(c(global_model$slab),seq(length(global_model$slab)))
names(slab_funnel)  <- c("name","id")

#reset all design settings
dev.off()

par(font=1,cex=0.8,srt=25)
funnel(global_model, xlim =c(-1.8,4.8),xlab = "Observed Outcome (Hedges'g)", ylab ="Standard Error (SE)",label = "all",slab=slab_funnel$id,offset =0.8)
par(font=1,cex=0.7,srt=0)
rect(2.4,0.43,4.9,-0.005,density = NA, col = "#ffffff",border="black")
text(2.6,seq(0.01,0.41,0.025), slab_funnel$id,adj=0)
text(2.9,seq(0.01,0.41,0.025), slab_funnel$name,adj=0)

# dev.print(pdf, file=eval(paste("Diagnostics/funnel_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))
# #fail safe n
# print(fsn(yi=g, vi = v,data=data_set)
# print(fsn(yi=g, vi = v,data=data_set,type="Orwin",target =0.2))
# 
# #Rank Correlation Test
# print(ranktest(res_rob))
# 


