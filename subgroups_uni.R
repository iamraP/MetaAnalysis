
#dev.new()
par(xpd=NA)

############# SETTINGS #############
correction_mode = 1
remove_outliers <- TRUE


### a little helper function to add Q-test ####
uni_hetero <- function(text, x) {
  bquote(paste(.(text), " (Q(", .(x$k - x$p),") = ", .(formatC(x$QE, digits=2, format="f")),
									   
               ", p ", .(metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep=" ")),
               "I²", " = ", .(formatC(x$I2, digits=1, format="f")), "%, ",
               "τ²", " = ", .(formatC(x$tau2, digits=2, format="f")), ")"))}

mv_hetero <- function(text, x, effects) {
  sigma_strings <- c()
  for (i in 1:length(effects)){
    sigma_strings<- c(sigma_strings, paste(", σ²_", effects[i]," = ", (formatC(x$sigma2[i], digits=2, format="f")),sep=""))
  }                  
  sigma_string <- paste(sigma_strings,collapse="")
  p_value <- metafor:::.pval(x$QEp, digits=2, showeq=TRUE, sep=" ")
  Q_value <- formatC(x$QE, digits=2, format="f")
  
  paste(text,
        " (Q(", x$k - x$p,") = ", Q_value ,
							 
        ", p ", p_value,
        sigma_string,")")}



##Get Data ####


setwd("G:/Meine Ablage/PhD/FMT_MetaAnalysis")
data <- read.csv("export_effect.csv",sep =",",fileEncoding="UTF-8-BOM")

data$index <- 1:nrow(data)

effects <- list("mod", "bsl_pre_post", "bsl_mod")
correction <- c("chen", "con", "lib")


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
# add name tag for the plot
data_set$names <- paste(data_set$Author,as.character(data_set$Year),sep=", ")

# drop studies withoout data, i.e. those w/o effectsizes
data_set <- drop_na(data_set)

###### calculate the overall effect ##########
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

## calculate  different subgroup effects  #########
#todo add the other subgroups below!

# get values for the rows of the forest plot
n_grp_rows <- c()
# keep track which studies were not included in the model

subsets <- list(### Comparisons
  subset(data_set,effect_type =="mod"),
  subset(data_set,effect_type =="bsl_pre_post"),
  subset(data_set,effect_type =="bsl_mod"),    #incl. outlier
  #no of sessions
  subset(data_set,number.of.sessions == 1), #single session
  subset(data_set,number.of.sessions > 1), #multiple   #incl. outlier
  # Modulation direction
  subset(data_set,Direction == -1),              # Down
  subset(data_set,Direction == 1),             #Up   #incl. outlier
  # individual frequencies             #
  subset(data_set,IF == 1),
  subset(data_set,IF == 0),             # standard freq    #incl. outlier
  # number of frequencies,    
  subset(data_set,n_freq > 1),         # Muti Band
  subset(data_set,n_freq == 1))                  # Single Band #incl. outlier      
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
outlier_subsets <- list()
outlier_subsets_labels <-list()
n_grp_rows <- c()
sink(file = "inf_out_dia_uni.txt") # store evaluation in text file 
for (k in 1:length(subsets)){
  sub_set <-subsets[[k]]
  cat("______________________________________________________________________________________\n______________________________________________________________________________________\n                              ",subsets_label[k], "\n______________________________________________________________________________________\n")
  print(subsets_label[k])
  if (k<=3){
    n_grp_rows <- c(n_grp_rows,nrow(sub_set))
  }
  res <- rma(yi = as.numeric(g),
             vi=as.numeric(v),
             data=sub_set,
             slab =names)
  cook <-cooks.distance.rma.uni(res)
  hat <- hatvalues.rma.uni(res)
  res_non_robust <- res


  # Hat values
  #print(hat)
  
  #Students residuals  (with non-robust model)
  stud_res <- residuals.rma(res_non_robust,type="rstudent")
  
  assign(paste("res",correction[correction_mode],subsets_label[k],sep = "."), res)
  
  
  
  outliers <- names(stud_res[sqrt(stud_res^2)>1.96])
  influencers <- names(cook[sqrt(cook^2)>median(cook)+6*IQR(cook)])         
  
  #print(stud_res)
  print(paste("cutoff for outliers: 1.96"))
  print(paste("outliers: ", outliers))
  
  #print(cook)
  print(paste("cutoff for influencers: ", median(cook)+6*IQR(cook)))
  print(paste("influencers: ", influencers))
  if (any(outliers %in% influencers)){
    outlier_subsets <-  append(outlier_subsets,list(subset(sub_set,!(sub_set$names %in% outliers[outliers %in% influencers]))))
    outlier_subsets_labels <-  append(outlier_subsets_labels, subsets_label[k])
  }

  
  
}
sink(file = NULL)


# ############### SUBSET OUTLIER EVALUATION - INVESTIGATE BY HAND!! ############## 

# # estimate variance
profile(res_rob,sigma=1)
dev.print(pdf, file=eval(paste("Diagnostics/profile_studyID_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))
profile(res_rob,sigma=2)
dev.print(pdf, file=eval(paste("Diagnostics/profile_effectType_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))



######## same analysis w/o the influential outlier #####

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

# #### subgroup models wihtout the outlier ### 
# outlier_subsets <- list(subset(data_set_cleaned,effect_type =="bsl_mod"),#incl. outlier
#                     subset(data_set_cleaned,number.of.sessions > 1), #multiple   #incl. outlier
#                     subset(data_set_cleaned,Direction == 1),                #Up   #incl. outlier  
#                     subset(data_set_cleaned,IF == 0),              # standard freq    #incl. outlier
#                     subset(data_set_cleaned,n_freq == 1))             # Single Band #incl. outlier
# 
# outlier_subsets_labels <- c("bsl_mod",
#                        "multi_sess",
#                        "up",
#                        "bsl_mod",
#                        "fix_freq",
#                        "single_band")


for (k in 1:length(outlier_subsets)){
  sub_set <- outlier_subsets[[k]]
  print(sum(duplicated(sub_set$study.ID)))
  if (sum(duplicated(sub_set$study.ID)) > 0) {
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
    res <- rma.mv(yi = as.numeric(g),
                  V = V,
                  random = ~ 1 | study.ID,
                  data=sub_set,
                  slab =names)
    res <- robust(res,
                  cluster = study.ID,
                  adjust = TRUE,
                  clubSandwich = TRUE)
  } else {
    res <- rma(yi = as.numeric(g),
               vi=as.numeric(v),
               data=sub_set,
               slab =names)
  }
  
  assign(paste("res_cleaned",correction[correction_mode],outlier_subsets_labels[k],sep = "."), res)
  
}

####### Settings for forest plot

# calculate rows required in the forest plot for studies not inclueded in the model 
n_grp_rows <- c(n_grp_rows,nrow(data) - sum(n_grp_rows))

#x- positions for the labels
#x_pos <- c(-8.8,-7.6,-6.8,-5.6,-4.8,-3.6,-2.8,seq(7.74,9.9,0.63))

x_pos <- c(-7.7,-6.9,-6.4,-5.6,-5.1,-4.3,-3.8,seq(7.74,9.9,0.63))
# row distribution

            # space header          # space subgroup model    #space for rows
req_rows <- length(n_grp_rows)*2 + (length(n_grp_rows)-1)*2+ sum(n_grp_rows)

y_lim=c(-24, req_rows)

# row distribution
poly_rows<- c() # ploygons for subgroups wihtin the plot
header_rows <- c(req_rows-3) # for the headers
add_rows <- c() # for studies not included in the analyis
forest_rows <- c() # for inclued studies
stop <- req_rows-4  # first studie
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


####### relable data 
data_set[data_set=="Across"] <- "A"
data_set[data_set=="Within"] <- "W"   
data_set$Direction[data_set$Direction==-1] <- "↓"
data_set$Direction[data_set$Direction==1] <- "↑"
data_set$IF[data_set$IF==0] <- "F"
data_set$IF[data_set$IF==1] <- "I"
data_set$add_freq[data_set$add_freq=="alpha"] <- "α"
data_set$add_freq[data_set$add_freq=="beta"] <- "β"
data_set[data_set=="Y"] <- "++"
data_set[data_set=="mY"] <- "+"
data_set[data_set=="mN"] <- "-"
data_set[data_set=="N"] <- "--"

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


############ Forrest Plot ################## 

#get the robust model of the currently applied correction 
res_rob <- eval(parse(text=paste("res_rob",correction[correction_mode],sep = ".")))

#mar:  bottom, left, top and right
par(oma =c(0,0,0,0), mar= c(6,0,0,12),font=1,cex=0.8)

### get the weights and format them as will be used in the forest plot
waldi <-forest(res_rob,
               xlim = c(-16.2,6.9),
               addpred = TRUE,
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
               mlab= eval(mv_hetero("REM for All Studies",res_rob,c("s","e"))),
               psize=0.9,
               header="Author(s) and Year")

#mar outlier with an asterix 
# text(-14.8,16.1,"*",cex=2)
# text(2.68,16.1,"*",cex=2)


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
}

#w/o outlier
addpoly(res_rob_cleaned, row=-2.4, mlab=eval(mv_hetero("REM w/o influential outlier",res_rob_cleaned,c("s","t"))),addpred=TRUE,col="#a4a4a4")



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
    #colorer cleaned models properly
    if(k==2){col<-"#a4a4a4"
      description <- "REM w/o influential outliers for"
    } else {col <-"#000000"
      description <- "REM for"}
    # cat("______________________________________________________________________________________ \n______________________________________________________________________________________")
    # print(subsets_label[i])
    # print(subgroup)
    if (is.null(subgroup$random)){
      addpoly(subgroup, y_pos_subgroups, mlab=eval(uni_hetero(paste(description, subsets_label[i], sep=" "),subgroup)),addpred=TRUE)
    } else if (grepl("effect_type",subgroup$random)) {
      addpoly(subgroup, y_pos_subgroups, mlab=eval(mv_hetero(paste(description, subsets_label[i], sep=" "),subgroup,c("s", "e"))),addpred=TRUE)
    } else{
      addpoly(subgroup, y_pos_subgroups, mlab=eval(mv_hetero(paste(description, subsets_label[i], sep=" "),subgroup,c("s"))),addpred=TRUE,col=col)
    }
  }
}

#add xlab
# add percentage sign
# add description
text(c(-7.8,-6.6,-5.3,-4,8.5),33, pos=3,cex=.8, c("Direction",
                                                  "Subjects",
                                                  "Sessions",
                                                  "Frequency",
                                                  "Study DIAD - Risk of Bias"))

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

text(0,-30, "Observed Outcome (Hedge's g)",font=1)

#segments(6, -30, x1 = 6, y1 = 32,lty = "dotdash")

text(c(rep(6.8,4),rep(7.8,5)),seq(-1,-9),pos =4, c(
  "(A) Fit between Concept & Operations",
  "(B) Clarity of Causal Inference",
  "(C) Generalizability of Findings",
  "(D) Precision of Outcome Estimation",
  "Ratings:",
  "++ : \"Yes\"",
  "+  : \"Maybe Yes\"",
  "-  : \"Maybe No\"",
  "-- : \"No\""),
  cex=0.8,font=1)


# ToDo include into string creation functions, add parameter for inf outlier, select of subset of symbols, this might be a lot of work 

# mark influential studies
#text(-14.1,8.2,cex=0.8, c("†")) #Chen
#text(-7.2,5.5,cex=0.8, c("†"))  # bsl mod
#text(-6.9,-6,cex=0.8, c("†"))#single
#Wang
#text(-13.95,16.4,cex=1.4, c("*")) #Wang
text(-6.8,-2,cex=1.4, c("*")) # overall model
#text(-7.35,5.7,cex=1.4, c("*"))  # bsl mod
#text(-7.1,-8.6,cex=1.4, c("*")) #multi
text(-7.3,-11.2,cex=1.4, c("*")) #up
text(-6.8,-16.7,cex=1.4, c("*")) #fix


#dev.print(pdf, file=eval(paste("Diagnostics/Forest_",correction[correction_mode],"_clean",remove_outliers,".pdf",sep="")))




