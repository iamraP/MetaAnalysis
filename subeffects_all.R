

################## Load  Data
data <- read.csv("export_effect.csv",sep =",",fileEncoding="UTF-8-BOM")

correction <- c("chen", "con", "lib")
correction_mode = 2
remove_outliers <- FALSE
################ Select correction mode
if (remove_outliers){
  data <- subset(data,data$study.ID!="RJHQE9FU")
}

data_set <-  data %>% select(-contains(correction[-correction_mode]))


hedges_g <- paste("g",correction[correction_mode],sep="_")
variance <-paste("var_g",correction[correction_mode],sep="_")
data_set <- data_set %>% rename("g" = hedges_g, "v" = variance)

data_set <- data_set[c("Author",
                       "Year",
                       "grp1","grp1_size",
                       "grp2","grp2_size",
                       "t1","t2",
                       "study.ID",
                       "g",
                       "v",
                       "effect_type",
                       "Direction",
                       "Session",
                       "PaperID",
                       "number.of.sessions",
                       "A",
                       "B",
                       "C",
                       "D")]


data_set$names <- paste(data_set$Author,as.character(data_set$Year),sep=", ")



# drop studies without data
data_set <- drop_na(data_set)


#################  Run  sub analysis for each of the following



#### Direction 
# Up
sub_set <- subset(data_set,Direction == -1)
#Down 
sub_set <- subset(data_set,Direction == 1)

#### Within/Across Sessions
sub_set <- subset(data_set,Session=="Within")
sub_set <- subset(data_set,Session=="Across") #run with and without outliers

#### Comparisons
sub_set <- subset(data_set,effect_type =="mod")
sub_set <- subset(data_set,effect_type =="bsl_pre_post")
sub_set <- subset(data_set,effect_type =="bsl_mod") # run with and without outliers


# Single vs Mulltiple Sessions:
sub_set <- subset(data_set,number.of.sessions == 1)
sub_set <- subset(data_set,number.of.sessions > 1) # run with and without outliers


############# Start Subanalysis - use multilevel model, whenever two effects from the same study are inclueded

sub_set$names <-  paste(sub_set$Author,as.character(sub_set$Year),sep=", ")

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

res

profile(res)


########## study diad results:
options <- list("Y", "mY", "mN", "N")#
items <- list(data$A, data$B,data$C,data$D)

for (i in 1:4){
  print("________________________________________")
  print(i)
  for (k in 1:4){
    print(paste(options[k],length(items[i][items[[i]]==options[k]])))
  }
}

