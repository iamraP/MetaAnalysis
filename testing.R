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
                       "g",
                       "v",
                       "effect_type",
                       "Direction",
                       "Subjects",
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

# drop studies withoout data, i.e. those without effectsizes
data_set <- drop_na(data_set)

test <- rma(yi = as.numeric(g),
    vi=as.numeric(v),
    data=data_set,
    slab =names)



(res_rob$QE^2-(res_rob$k-1))/res_rob$QE
