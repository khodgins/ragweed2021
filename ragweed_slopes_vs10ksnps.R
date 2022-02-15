
#####
# START HERE
#####

# read tables into R
library(data.table)
library(dplyr)

setwd("~/ha22/ragweed2/baypass/")

# North America
na.df = as.data.frame(fread("NA-2-NULL10k.gt", header = F))
rownames(na.df) = scan("NA-2-NULL10k.snps", what = "list")

na.pops.lat = distinct(read.table("NA-baypass-samples-pops-env.txt")[, -1])[, c(1,2)][sort(rep(1:43, 2)), ]


#get dim
no<-dim(na.df)[2]
snp<-dim(na.df)[1]


#get locations of each allele count
a1<-seq(1,no,2)
a2<-seq(2,no,2)

####turn into a loop for each snp
co<-vector()
for(i in 1:snp) {
  #get the allele counts for each pop
  A1<-t(na.df[i,a1])
  A2<-t(na.df[i,a2])
  na.pops.lat$V3[a1]
  
  #make a dataframe with the allele counts for each pop followed by the lat
  df<-cbind.data.frame(A1,A2, na.pops.lat$V3[a1])
  rownames(df)<-na.pops.lat$V2[a1]
  colnames(df)<-c("a1","a2", "lat")
  
  #run glm
  glm<-glm(data = df, cbind(a1, a2) ~ lat, family=binomial)
  #save slope for lat
  co[i]<-glm$coefficients[2]
}

#write.table
write.table(co, file="na_slope10k.txt")



# read tables into R
library(data.table)
library(dplyr)

setwd("~/ha22/ragweed2/baypass/")

# Europe
eu.df = as.data.frame(fread("EU-2-NULL10k.gt", header = F))
rownames(eu.df) = scan("EU-2-NULL10k.snps", what = "list")

eu.pops.lat = distinct(read.table("EU-baypass-samples-pops-env.txt")[, -1])[, c(1,2)][sort(rep(1:31, 2)), ]

# your allele count tables are na.df and eu.df
# the pop names and corresponding lats are in na.pops.lat and eu.pops.lat
# (entries are doubled to match columns in the allele tables)

#get dim
no<-dim(eu.df)[2]
snp<-dim(eu.df)[1]


#get locations of each allele count
a1<-seq(1,no,2)
a2<-seq(2,no,2)

####turn into a loop for each snp
co<-vector()
for(i in 1:snp) {
  #get the allele counts for each pop
  A1<-t(eu.df[i,a1])
  A2<-t(eu.df[i,a2])
  eu.pops.lat$V3[a1]
  
  #make a dataframe with the allele counts for each pop followed by the lat
  df<-cbind.data.frame(A1,A2, eu.pops.lat$V3[a1])
  rownames(df)<-eu.pops.lat$V2[a1]
  colnames(df)<-c("a1","a2", "lat")
  
  #run glm
  glm<-glm(data = df, cbind(a1, a2) ~ lat, family=binomial)
  #save slope for lat
  co[i]<-glm$coefficients[2]
}

#write.table
write.table(co, file="eu_slope10k.txt")


##########slope test versus 10k non-genic SNPs
#from ~/ha22/ragweed2/baypass/
na_slope<-read.table("na_slope10k.txt")
eu_slope<-read.table("eu_slope10k.txt")
hist(na_slope$x)
hist(eu_slope$x)

#Scaf 21 - highly significant in NA, only marginally sig in EU
merge.full.range <- read.csv("scaf_21_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x>glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)



#Scaf 5 - not sig
merge.full.range <- read.csv("scaf_5_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)

#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x<glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x>glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)


#not sig in model so don't test
merge.full.range <- read.csv("scaf_2_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)

#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x>glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)


#sig in NA not EU
merge.full.range <- read.csv("scaf_31_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x<glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)

#not in model so  - don't test
merge.full.range <- read.csv("scaf_448_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)

#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x>glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)

#27b sig in EU & NA
merge.full.range <- read.csv("scaf_27_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x>glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)

#27a sig in EU not NA
merge.full.range <- read.csv("scaf_27ii_map_model.csv", header = T)
merge.full.range$time<-as.factor(merge.full.range$time)
merge.full.range$range<-as.factor(merge.full.range$range)
#model with just lat (unadjusted)
glm<-glm(data = merge.full.range[merge.full.range$range=="na" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(na_slope$x[na_slope$x>glm$coefficients[2]])/length(na_slope$x)
#two sided
length(na_slope$x[abs(na_slope$x)>abs(glm$coefficients[2])])/length(na_slope$x)

glm<-glm(data = merge.full.range[merge.full.range$range=="eu" & merge.full.range$time=="modern",], cbind(h1count, h2count) ~ lat , family=binomial)
#one sided
#length(eu_slope$x[eu_slope$x<glm$coefficients[2]])/length(eu_slope$x)
#two sided
length(eu_slope$x[abs(eu_slope$x)>abs(glm$coefficients[2])])/length(eu_slope$x)

