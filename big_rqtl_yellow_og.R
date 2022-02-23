####QTL ANALYSIS ON BIG MAP - YELLOW FAMILY###
install.packages("qtl")
library(qtl)
setwd("~/Downloads/")
getwd() 


##load the csv file##
big_yellow_og_gp <- read.cross("csv", file = "big_yellow_rqtl_input.csv", genotypes = NULL, crosstype = "4way")


big_yellow_og_gp <- jittermap(big_yellow_og_gp)
summary(big_yellow_og_gp)
#4-way cross

#No. individuals:    336 

#No. phenotypes:     16 
#Percent phenotyped: 100 100 100 99.4 97.6 98.2 99.4 97 97 97.6 94.6 98.2 98.8 
#99.1 99.1 99.1 

#No. chromosomes:    18 
#Autosomes:      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 

#Total markers:      4493 
#No. markers:        513 489 488 164 365 186 328 65 312 307 301 106 219 195 155 
#139 63 98 
#Percent genotyped:  100 
#Genotypes (%):      AC:25.3  BC:26.6  AD:24.9  BD:23.2  AC/AD:0.0  BC/BD:0.0 
#AC/BC:0.0  AD/BD:0.0  AC/BD:0.0  BC/AD:0.0  not AC:0.0 
#not BC:0.0  not AD:0.0  not BD:0.0

##after removing 100 linked markers ~4393


##plot genetic map##
plotMap(big_yellow_og_gp)

##plot phenotypes to assess normality 
par(mfrow=c(3,6))
for(i in 10:14)
  plotPheno(big_yellow_og_gp,pheno.col=i)

par(mfrow=c(2,3))
for(i in 10:14)
  plotPheno(yellow_gp,pheno.col=i)
par(mfrow=c(2,3))
plotPheno(big_yellow_og_gp,pheno.col=10, xlab="Initial budding time (days)")
plotPheno(big_yellow_og_gp,pheno.col=11,xlab="1 cm budding time (days)")
plotPheno(big_yellow_og_gp,pheno.col=12,xlab="Male flowering time (days)")
plotPheno(big_yellow_og_gp,pheno.col=13,xlab="Female flowering time (days)")
plotPheno(big_yellow_og_gp,pheno.col=14,xlab="Height (cm)")

big_yellow_og_untransformed<- read.cross("csv", file = "big_yellow_rqtl_input.csv", genotypes = NULL, crosstype = "4way")
big_yellow_og_untransformed <- jittermap(big_yellow_og_untransformed)



par(mfrow=c(2,3))
for(i in 10:14)
  plotPheno(yellow_gp,pheno.col=i)
par(mfrow=c(2,3))
plotPheno(big_yellow_og_untransformed,pheno.col=10, main="a", ylab="Frequency", xlab="Initial budding time (days)")
plotPheno(big_yellow_og_untransformed,pheno.col=11, main="b", ylab="Frequency", xlab="1 cm budding time (days)")
plotPheno(big_yellow_og_untransformed,pheno.col=12, main="c", ylab="Frequency", xlab="Male flowering time (days)")
plotPheno(big_yellow_og_untransformed,pheno.col=13, main="d", ylab="Frequency", xlab="Female flowering time (days)")
plotPheno(big_yellow_og_untransformed,pheno.col=14, main="e", ylab="Frequency", xlab="Height (cm)")

#####log transform all but height
big_yellow_og_gp <- transformPheno(big_yellow_og_gp,pheno.col=(c(10:13)),transf=log)
par(mfrow=c(3,6))
for(i in 10:14)
  plotPheno(big_yellow_og_gp,pheno.col=i)


##SEGREGATION DISTORTION####
#use geno.table to check segregation against Mendelian proportions
#small p-value means segregation distortion - v different to expected under Mendelian
gt_filt_big_yellow_og<-geno.table(big_yellow_og_gp)
seg_dist_filt_big_yellow_og<-gt_filt_big_yellow_og[gt_filt_big_yellow_og$P.value<1e-4,]
View(seg_dist_filt_big_yellow_og)
str(seg_dist_filt_big_yellow_og)

#untransformed 
gt_filt_big_yellow_og_ut<-geno.table(big_yellow_og_untransformed)
seg_dist_filt_big_yellow_og_ut<-gt_filt_big_yellow_og_ut[gt_filt_big_yellow_og_ut$P.value<1e-7,]
View(seg_dist_filt_big_yellow_og_ut)
str(seg_dist_filt_big_yellow_og_ut)


#DROP markers showing SEG DIST ####
#make the markers id.d with high segregation distortion into an object
SD_markers_big_yellow_og <- row.names(seg_dist_filt_big_yellow_og)
str(SD_markers_big_yellow_og)

#untransformed 
SD_markers_big_yellow_og_ut <- row.names(seg_dist_filt_big_yellow_og_ut)
str(SD_markers_big_yellow_og_ut)

#drop these markers from the cross object
big_yellow_og_gp <- drop.markers(big_yellow_og_gp,SD_markers_big_yellow_og)
summary(big_yellow_og_gp)

#No. chromosomes:    17 
#Autosomes:      1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 

#Total markers:      4095
#Total markers:      3995 

plotMap(big_yellow_og_gp)

#untransformed 
big_yellow_og_untransformed <- drop.markers(big_yellow_og_untransformed,SD_markers_big_yellow_og_ut)
summary(big_yellow_og_untransformed)

##COMPARE INDIVIDUALS' GENOTYPES####
#are there any pairs of genotypes that are unusually similar? (mix-ups)
par(mfrow=c(1,1))
cg <- comparegeno(big_yellow_og_gp)
hist(cg, breaks=200,xlab="Proportion of identical genotypes")
rug(cg)
which(cg > 0.9, arr.ind=TRUE)
#none above 80 or 90 so okay 

save(big_yellow_og_gp,file="big_yellow_og_gp.RData")
load("big_yellow_og_gp.RData")


##CHECK FOR LINKAGE BETWEEN MARKERS ON DIFFERENT CHROMOSOMES####
#use est.rf to estimate the recombination fractions between all pairs of markers
#it inserts info back into the cross object
big_yellow_og_gp <-est.rf(big_yellow_og_gp)

#look at pairwise recombination fractions
plotRF(big_yellow_og_gp)
plotRF(big_yellow_og_gp, chr=c(14,15))
plotRF(big_yellow_og_gp, chr=c(2))

plotMap(big_yellow_og_gp)
#some markers on LG14 (scaffold 5) and 15 (scaffold 1) appear to be linked 

##check with these markers removed - all good 

save(big_yellow_og_gp,file="big_yellow_og_gp.RData")
load("big_yellow_og_gp.RData")
summary(big_yellow_og_gp)

##use RIPPLE to check marker order within chromosomes####
#recommended to use crossover count method with a large window
#then max likelihood with a smaller window
#for all chromosomes:
rip <- vector("list", nchr(big_yellow_og_gp))
names(rip) <- names(big_yellow_og_gp$geno)
for(i in names(big_yellow_og_gp$geno)) rip[[i]] <- ripple(big_yellow_og_gp, i, 7, verbose=FALSE)
#then work out the difference in initial crossovers vs best of the other orders
dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])
dif.nxo
#all came out as 0 or -1, indicating that the next best order is the same as this current order (or worse)
#this code below checks ALL alternate orders
any(dif.nxo > 0)
#FALSE --> there are no better orders

#try using the crossover method with a smaller window (3 then 2):
for(i in names(big_yellow_og_gp$geno)) rip[[i]] <- ripple(big_yellow_og_gp, i, 3, verbose=FALSE)
dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])
any(dif.nxo > 0)
#FALSE --> there are no better orders


##IDENTIFYING GENOTYPING ERRORS####
#search for unusually tight double cross-overs using genotyping error LOD scores to measure the evidence for genotyping errors
#should use a map based on the data, which for this data is already the case
#use calc.errorlod to calculate the gtype err LOD
big_yellow_og_gp_err <- calc.errorlod(big_yellow_og_gp)
#use top.errorlod to find genotypes with error LOD score higher than specified cut off 
#- use 5 as per guide
top <- top.errorlod(big_yellow_og_gp_err, cutoff=4)
#-->  No errorlods above cutoff.
#none above 5 or 3 either 

###COUNTING CROSSOVERS###
#use countXO to count the number of crossovers within each individual
#low or high = suspicious
nxo <- countXO(big_yellow_og_gp)
plot(nxo, ylab="No. crossovers")
nxo[nxo<3]
#no major outliers

##SAVE PROGRESS
save(big_yellow_og_gp,file="big_yellow_og_gp.RData")
load("big_yellow_og_gp.RData")

##RUN SINGLE-QTL SCANS####
#calculate genotype probabilities
big_yellow_og_gp <- calc.genoprob(big_yellow_og_gp,error.prob = 0.01)

#run the scan for BUDDING, FLOWERING and HEIGHT 
#use extended H-K method
big_yellow_og_gp.scan.out.ehk <- scanone(big_yellow_og_gp,pheno.col=c(10:14),method="ehk")
#5: In checkcovar(cross, pheno.col, addcovar, intcovar, perm.strata,  :
#Dropping 3-18 individuals with missing phenotypes.

save(big_yellow_og_gp.scan.out.ehk,file="big_yellow_og_gp.scan.out.ehk.RData")
load("big_yellow_og_gp.scan.out.ehk.RData")
par(mfrow=c(3,2))
for(i in 1:5)
  plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = i)
summary(big_yellow_og_gp.scan.out.ehk,threshold=3,format="allpheno")

#chr   pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#24962   1 21.65     3.69        2.668      4.642      2.864     2.007
#25673   1 56.16     3.17        2.504      3.468      3.900     3.246
#26289   1 79.92     1.52        1.498      1.262      1.245     5.710 
#9207    2 15.86    13.71       16.516     14.566     15.456     8.645 
#10159   2 93.17    18.58       19.808     18.735     20.347     6.687
#10167   2 93.40    18.45       21.081     18.507     20.349     6.994 
#10178   2 94.30    18.45       21.081     18.507     20.349     6.994
#18146   3 20.57     0.99        0.916      0.437      0.806     3.818
#7       4  3.52     2.03        3.138      2.096      1.908     0.394
#21481   6 37.97     2.59        2.446      3.668      4.361     1.164
#21543   6 40.92     2.71        2.757      3.661      4.523     1.130
#21842   6 55.65     1.09        1.677      1.811      2.793     3.238
#15616  12 23.49     5.11        4.736      5.180      6.637     0.318 


#The function scanone may also be used to perform a permutation test to get a genome-wide LOD significance threshold
big_yellow_og_gp.perm.out.ehk <- scanone(big_yellow_og_gp,pheno.col=c(10:14),method="ehk",n.perm=1000)
save(big_yellow_og_gp.perm.out.ehk,file="big_yellow_og_gp.perm.out.ehk.RData")
load("big_yellow_og_gp.perm.out.ehk.RData")
summary(big_yellow_og_gp.perm.out.ehk,alpha=c(0.001,0.05,0.1,0.2,0.5))
#LOD thresholds (1000 permutations)
#BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#0.1%     6.07         6.52       5.87       5.96      6.18
#5%       4.31         4.33       4.24       4.22      4.45
#10%      3.94         3.98       3.89       3.86      4.03
#20%      3.49         3.54       3.48       3.52      3.57
#50%      2.90         2.91       2.90       2.97      2.94

summary(big_yellow_og_gp.scan.out.ehk, perms=big_yellow_og_gp.perm.out.ehk,alpha=0.02,pvalues=TRUE)

      #chr pos BUD.DAYS pval BUD.1CM.DAYS  pval MFLOW.DAYS  pval FFLOW.DAYS  pval HEIGHT.CM  pval
#10159   2 93.2    18.58 0.00        19.81 0.000      18.74 0.000      20.35 0.000     6.687 0.001
#15616  12 23.5     5.11 0.01         4.74 0.025       5.18 0.006       6.64 0.001     0.318 1.000

summary(big_yellow_og_gp.scan.out.ehk, perms=big_yellow_og_gp.perm.out.ehk,alpha=0.001,pvalues=TRUE)
#only LG2 

####CREATE LOD PLOTS WITH SIGNIFICANCE THRESHOLDS (@ 0.1 and 5%)####
par(mfrow=c(2,3))
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 1, main="BUD.DAYS", ylab="LOD score", xlab="Linkage Group", alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.31,b=0,lty=2)
abline(a=6.07,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 2, main="BUD.1CM.DAYS", ylab="LOD score",alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.33,b=0,lty=2)
abline(a=6.52,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 3, main="MFLOW.DAYS", ylab="LOD score",alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.24,b=0,lty=2)
abline(a=5.87,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 4, main="FFLOW.DAYS", ylab="LOD score",alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.22,b=0,lty=2)
abline(a=5.96,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 5, main="HEIGHT.CM", ylab="LOD score",alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.45,b=0,lty=2)
abline(a=6.18,b=0,lty=2,col="gray")

#plot with a,b c instead of titles 
par(mfrow=c(2,3))
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 1, main="a", xlab="Linkage Group", ylab="LOD score",alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.31,b=0,lty=2)
abline(a=6.07,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 2, main="b", ylab="LOD score",xlab="Linkage Group", alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.33,b=0,lty=2)
abline(a=6.52,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 3, main="c", ylab="LOD score",xlab="Linkage Group", alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.24,b=0,lty=2)
abline(a=5.87,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 4, main="d", ylab="LOD score",xlab="Linkage Group", alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.22,b=0,lty=2)
abline(a=5.96,b=0,lty=2,col="gray")
plot(big_yellow_og_gp.scan.out.ehk,lodcolumn = 5, main="e", ylab="LOD score",xlab="Linkage Group", alternate.chrid=TRUE,ylim=c(0,25))
abline(a=4.45,b=0,lty=2)
abline(a=6.18,b=0,lty=2,col="gray")

#interval estimation of location of QTL#### 
#(Bayes credible intervals) - more consistent coverage so may be preferred

##BUD.DAYS (lodcolumn 1)
bayesint(big_yellow_og_gp.scan.out.ehk,chr=2,lodcolumn = 1, prob=0.95)
#chr      pos        BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#9930    2 74.39726 16.46245     18.43973   16.46561   18.45177  5.634301
#10159   2 93.16934 18.57597     19.80802   18.73544   20.34676  6.686804
#10240   2 94.97935 16.48883     19.10650   15.73207   18.33203  6.985039
bayesint(big_yellow_og_gp.scan.out.ehk,chr=12,lodcolumn=1,prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15665   12 25.30203 3.799745     3.731277   4.352865   4.980082 0.7392641

#BUD.1CM.DAYS
bayesint(big_yellow_og_gp.scan.out.ehk,chr=2,lodcolumn = 2, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#10019   2 81.13429 18.39781     20.57570   18.17627   20.28151  6.143304
#10178   2 94.30035 18.44614     21.08078   18.50713   20.34865  6.994279
#10240   2 94.97935 16.48883     19.10650   15.73207   18.33203  6.985039

bayesint(big_yellow_og_gp.scan.out.ehk,chr=12,lodcolumn = 2, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#15565  12 17.45801 2.826038     2.959147   2.854096   3.413077 0.5658858
#15616  12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15665  12 25.30203 3.799745     3.731277   4.352865   4.980082 0.7392641

#MFLOW.DAYS 
bayesint(big_yellow_og_gp.scan.out.ehk,chr=2,lodcolumn = 3, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#10019   2 81.13429 18.39781     20.57570   18.17627   20.28151  6.143304
#10159   2 93.16934 18.57597     19.80802   18.73544   20.34676  6.686804
#10178   2 94.30035 18.44614     21.08078   18.50713   20.34865  6.994279

bayesint(big_yellow_og_gp.scan.out.ehk,chr=12,lodcolumn = 3, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15665   12 25.30203 3.799745     3.731277   4.352865   4.980082 0.7392641

#FFLOW.DAYS 
bayesint(big_yellow_og_gp.scan.out.ehk,chr=2,lodcolumn = 4, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#9928    2 71.26026 16.45451     18.48334   16.48515   18.72458  5.423063
#10167   2 93.39534 18.44620     21.08078   18.50719   20.34869  6.994276
#10178   2 94.30035 18.44614     21.08078   18.50713   20.34865  6.994279

bayesint(big_yellow_og_gp.scan.out.ehk,chr=12,lodcolumn = 4, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15616   12 23.49002 5.108354     4.736440   5.179806   6.637249 0.3176855
#15636   12 24.28203 4.501926     4.519859   5.154904   5.978694 0.5787850

bayesint(big_yellow_og_gp.scan.out.ehk,chr=6,lodcolumn = 4, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#21312   6 24.69602 1.434972     1.620094   2.610052   2.919850  1.412087
#21543   6 40.92007 2.712353     2.756897   3.661052   4.523084  1.130354
#21967   6 62.67016 2.039872     2.747921   2.608480   3.549020  2.611866


#HEIGHT  *greater interval (similar to jlee)
bayesint(big_yellow_og_gp.scan.out.ehk,chr=2,lodcolumn = 5, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#9119    2 12.11601 13.79623     17.00363   14.82294   15.53793  7.638561
#9207    2 15.85602 13.71119     16.51626   14.56556   15.45572  8.644902 **
#10178   2 94.30035 18.44614     21.08078   18.50713   20.34865  6.994279

bayesint(big_yellow_og_gp.scan.out.ehk,chr=1,lodcolumn = 5, prob=0.95)
#chr      pos BUD.DAYS BUD.1CM.DAYS MFLOW.DAYS FFLOW.DAYS HEIGHT.CM
#26260   1 79.12945 1.561328     1.531650   1.272514   1.237558  5.155858
#26289   1 79.92145 1.521670     1.497589   1.261764   1.244562  5.709933
#26456   1 92.56651 1.358890     1.234544   1.106045   1.139035  4.339046

save(big_yellow_og_gp,file="big_yellow_og_gp.RData")
load("big_yellow_og_gp.RData")

###QTL EFFECT PLOTS####
par(mfrow=c(2,3))
for(i in 10:14)
  effectplot(big_yellow_og_gp,pheno.col = i, mname1="10159")
par(mfrow=c(2,3))
for(i in 10:15)
  effectplot(yell_filt,pheno.col = i, mname1="SC2ABSE_405.3_629873")
par(mfrow=c(2,3))
for(i in 10:15)
  effectplot(yell_filt,pheno.col = i, mname1="SC2ABSE_61.2_174096")

big_yellow_og_untransformed <- read.cross("csv", file = "/Users/dprapas/Desktop/QTL_inputs/BigMap_og/big_yellow_rqtl_input.csv", genotypes = NULL, crosstype = "4way")
big_yellow_og_untransformed <- jittermap(big_yellow_og_untransformed)

##EFFECT PLOTS 
#effect plots only for significant LOD traits & correct labels
par(mfrow=c(2,4))
effectplot(big_yellow_og_untransformed,pheno.col = 10, mname1="10159",xlab="Genotype",ylab="Initial budding time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 11, mname1="10159",xlab="Genotype",ylab="1 cm budding time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 12, mname1="10159",xlab="Genotype",ylab="Male flowering time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 13, mname1="10159",xlab="Genotype",ylab="Female flowering time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 14, mname1="10159",xlab="Genotype",ylab="Height (cm)")

par(mfrow=c(2,4))
effectplot(big_yellow_og_untransformed,pheno.col = 10, mname1="15616",xlab="Genotype",ylab="Initial budding time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 11, mname1="15616",xlab="Genotype",ylab="1 cm budding time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 12, mname1="15616",xlab="Genotype",ylab="Male flowering time (days)")
effectplot(big_yellow_og_untransformed,pheno.col = 13, mname1="15616",xlab="Genotype",ylab="Female flowering time (days)")


##PXG plots
par(mfrow=c(2,3))
for(i in 10:14)
  plotPXG(big_yellow_og_untransformed,pheno.col = i, marker="10159")

par(mfrow=c(2,3))
for(i in 10:13)
  plotPXG(big_yellow_og_untransformed,pheno.col = i, marker="15616")

par(mfrow=c(2,3))
for(i in 13)
  plotPXG(big_yellow_og_untransformed,pheno.col = i, marker="26289")


### RUN 2D, two-QTL scans#### This allows us to identify potential interactions among QTL and to assess evidence for linked QTL
## ^ takes a very long time, need to do via HPC (see scantwo_r_M3)
##upload combined scantwo and permutation test results

##FLOWERING traits 
load("big_yellow_og.FL.H.out2.RData")
load("big_yellow_og.FL.H.perm.all.hk.RData")

summary(big_yellow_og.FL.H.perm.all.hk)
#MFLOW.DAYS (1000 permutations)
#full   fv1  int  add  av1  one
#5%  12.7 10.20 8.78 7.46 4.73 4.15
#10% 12.2  9.78 8.34 7.11 4.43 3.87

#FFLOW.DAYS (1000 permutations)
#full  fv1  int  add  av1  one
#5%  12.6 10.2 8.87 7.62 4.85 4.33
#10% 12.0  9.7 8.34 7.12 4.43 3.93

#HEIGHT.CM (1000 permutations)
#full   fv1  int  add  av1  one
#5%  12.5 10.30 8.78 7.79 4.93 4.22
#10% 12.2  9.82 8.33 7.25 4.53 3.89

#MFLOW
summary(big_yellow_og.FL.H.out2, perms=big_yellow_og.FL.H.perm.all.hk, alphas=0.2,pvalues=TRUE,lodcolumn = 1)
##for 20%
         #pos1f pos2f lod.full  pval lod.fv1  pval lod.int pval   pos1a pos2a lod.add  pval lod.av1  pval
#c1:c2   21.6  93.4     26.3 0.000    7.56 0.976    3.03    1      21.6  93.2   23.27 0.000    4.53 0.081
#c1:c6   21.6  43.1     12.2 0.090    7.60 0.966    3.43    1      21.6  38.0    8.81 0.003    4.17 0.162
#c1:c12  21.6  23.5     11.5 0.244    6.29 1.000    1.78    1      21.6  23.5    9.69 0.002    4.51 0.085
#c2:c12  93.2 23.5     28.6 0.000    9.88 0.084    3.93    1      93.2  23.5   24.68 0.000    5.94 0.006
##for 5%
summary(big_yellow_og.FL.H.out2, perms=big_yellow_og.FL.H.perm.all.hk, alphas=0.05,pvalues=TRUE,lodcolumn = 1)
         #pos1f pos2f lod.full pval lod.fv1  pval lod.int pval     pos1a pos2a lod.add pval lod.av1  pval
#c2:c12  93.2  23.5     28.6    0    9.88 0.084    3.93    1      93.2  23.5    24.7    0    5.94 0.006

##FFLOW 
summary(big_yellow_og.FL.H.out2, perms=big_yellow_og.FL.H.perm.all.hk, alphas=0.2,pvalues=TRUE,lodcolumn = 2)
##for 20%
       #pos1f pos2f lod.full pval lod.fv1  pval lod.int pval     pos1a pos2a lod.add pval        lod.av1  pval
#c1:c2   5.23  90.5     28.8    0    8.40 0.584    3.58    1      55.9  93.4    25.2    0    4.82 0.053
#c2:c6  83.40  39.2     29.6    0    9.26 0.211    4.12    1      93.2  44.4    25.5    0    5.14 0.023
#c2:c12 93.2 23.5        32.4    0   12.07 0.001    3.65    1      93.2  23.5    28.8    0    8.42 0.000

summary(big_yellow_og.FL.H.out2, perms=big_yellow_og.FL.H.perm.all.hk, alphas=0.05,pvalues=TRUE,lodcolumn = 2)
#for 5%
#pos1f pos2f lod.full pval lod.fv1  pval lod.int pval             pos1a pos2a lod.add pval lod.av1  pval
#c2:c6   83.4  39.2     29.6    0    9.26 0.211    4.12    1      93.2  44.4    25.5    0    5.14 0.023 ********
#c2:c12  93.2  23.5     32.4    0   12.07 0.001    3.65    1      93.2  23.5    28.8    0    8.42 0.000

##HEIGHT 
summary(big_yellow_og.FL.H.out2, perms=big_yellow_og.FL.H.perm.all.hk, alphas=0.2,pvalues=TRUE,lodcolumn = 3)
##for 20%
#pos1f pos2f lod.full  pval lod.fv1  pval lod.int  pval            pos1a pos2a lod.add  pval lod.av1  pval
#c2:c6   15.9  56.1     14.8 0.004    6.15 1.000    1.40 1.000      15.9  56.1   13.39 0.000   4.744 0.067
#c2:c18  98.6   0.0     17.2 0.000    8.51 0.563    7.94 0.189      16.1  13.8    9.21 0.005   0.568 1.000

par(mfrow=c(3,2))
for(i in 1:3)
  plot(big_yellow_og.FL.H.out2,lodcolumn = i)
summary(big_yellow_og.FL.H.out2,threshold=3,format="allpheno")

plot(big_yellow_og.FL.H.out2, chr=c(1,2), lower = "cond-int")

##BUDDING traits 
load("big_yellow_og.BUD.out2.RData")
load("big_yellow_og.BUD.perm.all.hk.RData")


plot(big_yellow_og.BUD.out2, chr=c(2,12), lower = "cond-int")

plot(big_yellow_og.BUD.out2, chr=2, lower="cond-int", upper="cond-add")

summary(big_yellow_og.BUD.perm.all.hk)
#BUD.DAYS (1000 permutations)
       #full   fv1    int   add   av1   one
#5%     12.7   10.32  8.95  7.60  4.64  4.27
#10%    12.2   9.84   8.41  7.03  4.26  3.94

#BUD.1CM.DAYS (1000 permutations)
      #full   fv1   int   add   av1   one
#5%    12.5  10.09  8.70  7.43  4.74  4.23
#10%   12.0  9.68   8.37  7.05  4.40  3.86

#Tf,Tfv1(conditional-interactive),Ti,Ta,Tav1(conditional-additive)
#Note that (pos1f, pos2f) indicate the estimated positions of the QTL under the full model
#while (pos1a, pos2a) are the estimated positions under the additive model

#BUD
summary(big_yellow_og.BUD.out2, perms=big_yellow_og.BUD.perm.all.hk, alphas=0.2,pvalues=TRUE,lodcolumn = 1)
##for 20%
         #pos1f pos2f lod.full pval lod.fv1  pval   lod.int pval       pos1a pos2a lod.add pval lod.av1  pval
#c2:c12  93.2  23.5     28.3    0    9.75   0.116     3.32    1         93.2 23.5      25    0    6.42 0.002
summary(big_yellow_og.BUD.out2, perms=big_yellow_og.BUD.perm.all.hk, alphas=0.05,pvalues=TRUE,lodcolumn = 1)
##for 5%
        #pos1f pos2f lod.full pval lod.fv1  pval     lod.int pval            pos1a pos2a lod.add pval lod.av1  pval
#c2:c12  93.2  23.5     28.3    0    9.75   0.116    3.32    1             93.2  23.5      25    0    6.42 0.002*
#no interaction, only additive 

#BUD1CM 
summary(big_yellow_og.BUD.out2, perms=big_yellow_og.BUD.perm.all.hk, alphas=0.2,pvalues=TRUE,lodcolumn = 2)
##for 20%       
        #pos1f pos2f lod.full pval lod.fv1  pval lod.int pval     pos1a pos2a lod.add pval lod.av1  pval
#c2:c2   81.4  94.3     27.6    0    6.54 1.000    2.21    1      12.5  81.8    25.4    0    4.33 0.114
#c2:c3   93.8  59.1     28.6    0    7.54 0.956    3.34    1      93.2  55.7    25.3    0    4.20 0.137
#c2:c6   81.6  64.4     29.5    0    8.45 0.546    3.75    1      93.2  77.1    25.8    0    4.70 0.054
#c2:c12  93.2  23.5     30.8    0    9.74 0.092    4.03    1      93.2  23.5    26.8    0    5.71 0.006 
summary(big_yellow_og.BUD.out2, perms=big_yellow_og.BUD.perm.all.hk, alphas=0.05,pvalues=TRUE,lodcolumn = 2)
##for 5% 
         #pos1f pos2f lod.full pval lod.fv1  pval lod.int pval     pos1a pos2a lod.add pval lod.av1  pval
#c2:c12  81.6  23.5     30.8    0    9.74 0.092    4.03    1       93.2  23.5    26.8    0    5.71 0.006 *


find.marker(big_yellow_og_gp, 6, 44.4)
#[1] "21648"


#####MULTI-QTL FIT ####
big_yellow_og_gp<-sim.geno(big_yellow_og_gp,n.draws=1000,error.prob = 0.01)

big_yellow_og_untransformed<-sim.geno(big_yellow_og_untransformed,n.draws=1000,error.prob = 0.01)
##try H-K method (see rqtl notes) and back transform the log transformed phenotypes 

##FLOWERING

#FFLOW
#scanone results 
qtl.ff1 <- makeqtl(big_yellow_og_gp, chr=c(2, 12), pos=c(93.2,23.5))
big_yellow_og_FFLOW_out.fq <- fitqtl(big_yellow_og_gp, pheno.col=13,qtl=qtl.ff1, formula=y~Q1+Q2,get.ests=TRUE)
summary(big_yellow_og_FFLOW_out.fq)
plot(qtl.ff1)

qtl.ff1.ut <- makeqtl(big_yellow_og_untransformed, chr=c(2, 12), pos=c(93.2,23.5))
big_yellow_og_FFLOW_out.fq.ut <- fitqtl(big_yellow_og_untransformed, pheno.col=13,qtl=qtl.ff1.ut, formula=y~Q1+Q2,get.ests=TRUE)
summary(big_yellow_og_FFLOW_out.fq.ut)
plot(qtl.ff1)

#fitqtl summary

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 332 

#Full model result
----------------------------------  
#  Model formula: y ~ Q1 + Q2 

#df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
#Model   6  4.614223 0.76903721 27.98123 32.16741            0         0
#Error 325  9.730180 0.02993901                                         
#Total 331 14.344403                                                    


#Drop one QTL at a time ANOVA table: 
  ----------------------------------  
     #  df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
#2@93.2   3       3.352 21.343 23.370   37.32            0   < 2e-16 ***
#12@23.5  3       1.088  7.643  7.586   12.12            0  1.55e-07 ***
  ---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Estimated effects:
  -----------------
#  est       SE       t
#Intercept   4.11571  0.02472 166.502
#2@93.2.BC   0.02348  0.02752   0.853
#2@93.2.AD  -0.18427  0.02555  -7.211
#2@93.2.BD  -0.19844  0.02694  -7.367
#12@23.5.BC -0.00301  0.02600  -0.116
#12@23.5.AD  0.10190  0.02818   3.616
#12@23.5.BD  0.12453  0.02805   4.440

#scantwo results
qtl.ff2 <- makeqtl(big_yellow_og_gp, chr=c(2, 12, 6), pos=c(93.2,23.5,44.4))

rqtl.ff <- refineqtl(big_yellow_gp, qtl=qtl.ff2, formula=y~Q1+Q2+Q3,verbose=FALSE)
bayesint(rqtl.ff, qtl.index = 3)
#chr      pos      lod
#21199   6 14.99500 3.170072
#21200   6 16.16001 3.172300
#21992   6 70.86216 1.377265

big_yellow_og_FFLOW_out.fq2 <- fitqtl(big_yellow_og_gp, pheno.col=13,qtl=qtl.ff2, formula=y~Q1+Q2+Q3,get.ests=TRUE)
summary(big_yellow_og_FFLOW_out.fq2)

#fitqtl summary

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 332 

#Full model result
----------------------------------  
#  Model formula: y ~ Q1 + Q2 + Q3 

#df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
#Model   9  5.005387 0.55615415 30.93931 34.89436            0         0
#Error 322  9.339016 0.02900315                                         
#Total 331 14.344403                                                    


#Drop one QTL at a time ANOVA table: 
  ----------------------------------  
      #  df  Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
#2@93.2   3      3.3228 21.944 23.164  38.189        0.000   < 2e-16 ***
#12@23.5  3      0.7356  5.466  5.128   8.454        0.000  2.01e-05 ***
#6@44.4   3      0.3912  2.958  2.727   4.496        0.003   0.00416 ** 
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Estimated effects:
  -----------------
#  est       SE       t
#Intercept   4.18365  0.03171 131.924
#2@93.2.BC   0.01251  0.02729   0.458
#2@93.2.AD  -0.18908  0.02530  -7.475
#2@93.2.BD  -0.20283  0.02658  -7.632
#12@23.5.BC -0.00765  0.02639  -0.290
#12@23.5.AD  0.07938  0.02842   2.793
#12@23.5.BD  0.10515  0.02858   3.679
#6@44.4.BC  -0.03462  0.02699  -1.283
#6@44.4.AD  -0.08826  0.02747  -3.212
#6@44.4.BD  -0.08429  0.02911  -2.895

##MFLOW
qtl.mf1 <- makeqtl(big_yellow_og_gp, chr=c(2, 12), pos=c(93.2,23.5))
big_yellow_og_MFLOW_out.fq <- fitqtl(big_yellow_og_gp, pheno.col=12,qtl=qtl.mf1, formula=y~Q1+Q2,get.ests=TRUE)
summary(big_yellow_og_MFLOW_out.fq)

#fitqtl summary

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 330 

#Full model result
----------------------------------  
#  Model formula: y ~ Q1 + Q2 

#df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
#Model   6  2.919371 0.48656187 24.58075 29.03802            0         0
#Error 323  7.134246 0.02208745                                         
#Total 329 10.053617                                                    


#Drop one QTL at a time ANOVA table: 
  ----------------------------------  
#  df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
#2@93.2   3      2.2185 19.403 22.067  33.481            0   < 2e-16 ***
#12@23.5  3      0.6075  5.856  6.043   9.168            0  7.73e-06 ***
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Estimated effects:
  -----------------
#  est        SE       t
#Intercept   4.185626  0.021267 196.814
#2@93.2.BC   0.020738  0.023657   0.877
#2@93.2.AD  -0.149396  0.022080  -6.766
#2@93.2.BD  -0.161271  0.023195  -6.953
#12@23.5.BC -0.001022  0.022380  -0.046
#12@23.5.AD  0.079792  0.024203   3.297
#12@23.5.BD  0.092029  0.024167   3.808


#HEIGHT
qtl.h1 <- makeqtl(big_yellow_og_gp, chr=c(2), pos=c(93.2))
big_yellow_og_H_out.fq <- fitqtl(big_yellow_og_gp, pheno.col=14,qtl=qtl.h1, formula=y~Q1,get.ests=TRUE)
summary(big_yellow_og_H_out.fq)

#fitqtl summary (for 93.2)

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 333 

#Full model result
----------------------------------  
#  Model formula: y ~ Q1 

#df        SS        MS      LOD     %var Pvalue(Chi2)    Pvalue(F)
#Model   3  531.5682 177.18941 6.682451 8.827228 9.485521e-07 1.109754e-06
#Error 329 5490.3477  16.68799                                            
#Total 332 6021.9159                                                      


#Estimated effects:
  -----------------
#  est      SE      t
#Intercept 13.2211  0.4237 31.203
#2@93.2.BC  0.9943  0.6413  1.550
#2@93.2.AD -2.2001  0.5992 -3.672
#2@93.2.BD -1.6809  0.6341 -2.651


##BUDDING
qtl.bud <- makeqtl(big_yellow_og_gp, chr=c(2, 12), pos=c(93.2,23.5))

##BUD.DAYS 
big_yellow_og_BUD_out.fq <- fitqtl(big_yellow_og_gp, pheno.col=10,qtl=qtl.bud, formula=y~Q1+Q2,get.ests=TRUE)
summary(big_yellow_og_BUD_out.fq)

#fitqtl summary

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 328 

#Full model result
#----------------------------------  
#  Model formula: y ~ Q1 + Q2 

#df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
#Model   6  4.126976 0.68782937 24.43405 29.04028            0         0
#Error 321 10.084236 0.03141507                                         
#Total 327 14.211212                                                    


#Drop one QTL at a time ANOVA table: 
  ----------------------------------  
#  df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
#2@93.2   3      3.1433 19.325 22.118  33.352            0   < 2e-16 ***
#  12@23.5  3      0.8662  5.869  6.095   9.191            0  7.52e-06 ***
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Estimated effects:
  -----------------
#  est         SE       t
#Intercept   3.9678853  0.0254191 156.099
#2@93.2.BC   0.0279091  0.0284198   0.982
#2@93.2.AD  -0.1791197  0.0264556  -6.771
#2@93.2.BD  -0.1890806  0.0277445  -6.815
#12@23.5.BC -0.0003721  0.0266911  -0.014
#12@23.5.AD  0.0978117  0.0287896   3.397
#12@23.5.BD  0.1096363  0.0291279   3.764

##BUD.1CM
big_yellow_og_BUD.1CM_out.fq <- fitqtl(big_yellow_og_gp, pheno.col=11,qtl=qtl.bud, formula=y~Q1+Q2,get.ests=TRUE)
summary(big_yellow_og_BUD_out.fq)

#fitqtl summary

#Method: multiple imputation 
#Model:  normal phenotype
#Number of observations : 328 

#Full model result
----------------------------------  
#  Model formula: y ~ Q1 + Q2 

#df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
#Model   6  4.126976 0.68782937 24.43405 29.04028            0         0
#Error 321 10.084236 0.03141507                                         
#Total 327 14.211212                                                    


#Drop one QTL at a time ANOVA table: 
  ----------------------------------  
#  df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
#2@93.2   3      3.1433 19.325 22.118  33.352            0   < 2e-16 ***
# 12@23.5  3      0.8662  5.869  6.095   9.191            0  7.52e-06 ***
  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Estimated effects:
  -----------------
#  est         SE       t
#Intercept   3.9678853  0.0254191 156.099
#2@93.2.BC   0.0279091  0.0284198   0.982
#2@93.2.AD  -0.1791197  0.0264556  -6.771
#2@93.2.BD  -0.1890806  0.0277445  -6.815
#12@23.5.BC -0.0003721  0.0266911  -0.014
#12@23.5.AD  0.0978117  0.0287896   3.397
#12@23.5.BD  0.1096363  0.0291279   3.764


