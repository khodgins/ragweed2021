#is the change in slopes correlated with the modern estimates of the strength of selection
#read slope data - from sup mat table
slopes <- read.table("HB_slopes.txt", header = T)
#use full model
full<-slopes[slopes$model=="f",]
full$group<-ifelse(full$time=="historic" & full$range=="eu", 1,0)

#all HBs

#are hist eu slopes lower than ones that might be closer to equilibrium 
#should only use the ones with significant slopes since don't expect slope change if not under lat selection
#wilcox.test(abs(s) ~ group, data=full)

modna<-full[full$time=="modern" & full$range=="na",]
hisna<-full[full$time=="historic" & full$range=="na",]
modeu<-full[full$time=="modern" & full$range=="eu",]
hiseu<-full[full$time=="historic" & full$range=="eu",]

#modern na vs modern eu -> expect hist and mod to be correlated
cor.test(modna$s, modeu$s, method=c("pearson"))
cor.test(modna$slope, modeu$slope, method=c("pearson"))
plot( modna$slope, modeu$slope)
#change in slope in eu correlated with s in mod eu or mod na
cor.test((modeu$slope-hiseu$slope), modeu$s, method=c("pearson"))
plot(modeu$s, (modeu$slope-hiseu$slope), xlab="relative strength of selection", ylab="change in slope over time", col="red", pch=16)
cor.test((modna$slope-hisna$slope), modna$s, method=c("pearson"))
plot(modna$s, (modna$slope-hisna$slope), xlab="relative strength of selection", ylab="change in slope over time")


par(mfrow=c(2,2)) 
plot(modeu$s, (modeu$slope-hiseu$slope), xlab="relative strength of selection", ylab="change in slope over time",  pch=16, main="Europe")
plot(modna$s, (modna$slope-hisna$slope), xlab="relative strength of selection", ylab="change in slope over time",  pch=16, main="North America")

#only sig HBs
fulls<-full[full$lat_sig=="y",]
#are hist eu slopes lower than ones that might be closer to equilibrium
t.test(abs(slope) ~ group, data=fulls)
t.test(abs(fulls$slope[fulls$range=="eu"]) ~ fulls$time[fulls$range=="eu"])
t.test(abs(fulls$slope[fulls$range=="na"]) ~ fulls$time[fulls$range=="na"])

