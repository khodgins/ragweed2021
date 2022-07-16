#order forced figure 5 Diana's manuscript

library(ggplot2)
```{r}
x=34000000
y=2000000
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced27.mapped")
#male map
A<-ggplot(data=map, aes(x=V2, y=V3))+
  xlim(y, x)+ 
  ylim(0, 60)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 60,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 60,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height

#female map
B<-ggplot(data=map, aes(x=V2, y=V4))+
  xlim(y, x)+ 
  ylim(0, 40)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 40,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 40,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced27.mapped")
#male map
C<-ggplot(data=map, aes(x=V2, y=V3))+
  xlim(y, x)+ 
  ylim(0, 40)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 40,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 40,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height

#female map
D<-ggplot(data=map, aes(x=V2, y=V4))+
  xlim(y, x)+ 
  ylim(0, 45)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 45,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 45,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_27forced.mapped")
#male map
E<-ggplot(data=map, aes(x=V2, y=V3))+
  xlim(y, x)+ 
  ylim(0, 70)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 70,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 70,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height

#female map
F<-ggplot(data=map, aes(x=V2, y=V4))+
  xlim(y, x)+ 
  ylim(0, 45)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 45,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 45,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height
library(ggpubr)
ggarrange(A, B,C,D,E,F, 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3,common.legend = TRUE)
````
############
















#sex averaged
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/order27yellow2out.test2_gen_phy.txt")
plot(map$V2, map$V6)
abline(v=c(18930440,23019032), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))
abline(v=c(23966548,29469625), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))


#F1 unforced
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_448_2.map")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_2_2.map")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_31_2.map")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_21_2.map")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_27_2.map", sep="\t")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))




#########2
#forced
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced2.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_2forced.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced2.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink2.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


#########3
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced5.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_5forced.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced5.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink5.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


######21
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced21.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_21forced.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced21.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))



map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink21.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


#######31
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced31.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_31forced.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced31.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))



map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink31.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

####448

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced448.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_448forced.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))



map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced448.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink448.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))



###############27
#order not forced
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink27.mapped")
plot(map$V2, map$V3)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time
plot(map$V2, map$V4)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time
install.packages("scales")
library(scales)












#female map
ggplot(data=map, aes(x=V2, y=V4))+
  xlim(4000000, 31000000)+ 
  ylim(0, 50)+
  geom_point()+
  annotate("rect", xmin = 18930440, xmax = 23019032, ymin = 0, ymax = 50,
           alpha = .1,fill = "blue")+ #23019032
  annotate("rect", xmin = 23966548, xmax = 29469625, ymin = 0, ymax = 50,
           alpha = .1,fill = "purple")+  #hb27b
  labs(
    x = "chromosome position (bp)",
    y = "genetic distance (cM)",
  ) +
  theme_classic()+
  geom_vline(xintercept =  26795890, colour="black", size=2) + 
  geom_vline(xintercept=c(24111324,26860686), colour="grey", linetype=2) + #female flowering time
  geom_vline(xintercept=c(25529687,26860686), colour="grey", linetype=3) + ##male flowering, smallest interval
  geom_vline(xintercept=c(4229363,26860686), colour="grey", linetype=4)  #height




map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced27.mapped")
#male map
geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="black", alpha=0.5) +
plot(map$V2, map$V3, xlim=c(4000000, 31000000), ylim=c(0, 60), xlab="chromosome position (bp)", ylab="genetic distance (cM)"))
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA))
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height

#female map
plot(map$V2, map$V4, xlim=c(4000000, 31000000), ylim=c(0, 60), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced27.mapped")
#map map
plot(map$V2, map$V3, xlim=c(4000000, 31000000), ylim=c(0, 40), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height
#female map
plot(map$V2, map$V4, xlim=c(4000000, 31000000), ylim=c(0, 40), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_27forced.mapped")
#female
plot(map$V2, map$V3, xlim=c(4000000, 31000000), ylim=c(0, 80), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height
#male
plot(map$V2, map$V4, xlim=c(4000000, 31000000), ylim=c(0, 80), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mappinkforced27.mapped")
#female
plot(map$V2, map$V3,  xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height
#male
plot(map$V2, map$V4, xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/mapyellowforced27.mapped")
#female
plot(map$V2, map$V3,  xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height
#male
plot(map$V2, map$V4, xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_27forced.mapped")
#female
plot(map$V2, map$V3,  xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height
#male
plot(map$V2, map$V4, xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height



map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderF1_27_2.map", sep="\t")
plot(map$V2, map$V3, xlim=c(4000000, 31000000), ylim=c(0, 40), xlab="chromosome position (bp)", ylab="genetic distance (cM)")
rect(c(18930440), -1e6, c(23019032), 1e6, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1), border=NA)
rect(c(23966548), -1e6, c(29469625), 1e6, col=rgb(red = 0, green = 1, blue = 0, alpha = 0.1), border=NA)
abline(v=c(26795890), col=c("black"), lty=c(1), lwd=3) #central marker QTL-2
abline(v=c(24111324,26860686), col=c("grey","grey"), lty=c(2,2 ) , lwd=2)#female flowering time
abline(v=c(25529687,26860686), col=c("grey","grey"),lty=c(3,3 ), lwd=2)  #male flowering, smallest interval
abline(v=c(4229363,26860686), col=c("grey","grey"),lty=c(4,4 ), lwd=2) #height



ID<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/bigmap_phasematch.2.mapped2")
#from big_rqtl_yellow_og.R
big_yellow_og_gp$geno$`2`$map
pos<-big_yellow_og_gp$geno$`2`$map
mapped_diana<-ID[ ID$V3 %in% colnames(pos),]

plot(mapped_diana$V2, mapped_diana$V4, xlab="chromosome position (bp)", ylab="genetic distance (cM)")
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time, largest interval
abline(v=c(25529687,26795890,26860686), col=c("pink","pink","pink"))  #male flowering, smallest interval


##not helpful


##not helpful
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/F1_27.mapped")
plot(map$V2, map$V3)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time
plot(map$V2, map$V4)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time


#order not forced
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow27.mapped")

plot(map$V2, map$V3)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time
plot(map$V2, map$V4)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
abline(v=c(24111324,26795890,26860686), col=c("black","black","black")) #central marker female flowering time


#not forced
map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink27.mapped")
plot(map$V2, map$V3)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
plot(map$V2, map$V4)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))




map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink2.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink5.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink21.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink31.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderpink448.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow_nobad_27.mapped")
plot(map$V2, map$V3)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))
plot(map$V2, map$V4)
abline(v=c(18930440,23019032), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
abline(v=c(23966548,29469625), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))






map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow2.mapped")
plot(map$V2, map$V3)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(54075633,56509868), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow5.mapped")
plot(map$V2, map$V3)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(42452580,47751788), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow21.mapped")
plot(map$V2, map$V3)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(36078313,49457003), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))


map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow31.mapped")
plot(map$V2, map$V3)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(41057374,47742017), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))

map<-read.table("/Users/kayhodgins/Dropbox/Documents/ragweed/orderyellow448.mapped")
plot(map$V2, map$V3)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))
plot(map$V2, map$V4)
abline(v=c(17719474,33802585), col=c("blue", "blue"), lty=c(1,1), lwd=c(1, 1))



install.packages("MareyMap")
install.packages("tcltk")
install.packages("tkrplot")
install.packages("tools")
library(MareyMap)
startMareyMapGUI()

install.packages("onemap")
library(onemap)

data<-onemap_read_vcfR(
  vcf="/Users/kayhodgins/Dropbox/Documents/ragweed/yellow.fixed.grand.recode.vcf.gz",
  cross = "f2 intercross",
  parent1 = "G4",
  parent2 = "A3")

data.bins<-find_bins(data, exact=T)
create_data_bins(data, data.bins)
data.binned<-create_data_bins(data, data.bins)
test_segregation(data.binned)
data.seg_test<-test_segregation(data.binned)
plot(data.seg_test)
data.no_dist <- select_segreg(data.seg_test, distorted = FALSE, numbers = TRUE) 
data.twopts <- rf_2pts(data.binned, LOD = 4, max.rf = 0.4)

unique(bins_example$CHROM)
CHR1 <- make_seq(twopts, "ScBFxKa_27")
CHR1
CHR_mks <- group_seq(input.2pts = data, seqs = list(CHR1=ScBFxKa_27), 
                     unlink.mks = mark_no_dist, repeated = FALSE)

make_seq(data, "1")
data.LGs_LOD5<-group(data.mrks.no_dist, LOD = 5, max.rf = 0.5 )

set_map_fun(type="kosambi")
data.LG1 <- make_seq(data.LGs_LOD5, 1)
data.LG1.ord <- order_seq(data.LG1, n.init=5, touchdown=TRUE)
make_seq(data.LG1.ord, "force")
order_seq(data.LG1, n.init=5, touchdown=TRUE)

CHR2_ord <- order_seq(data$sequences)
CHR2_frame <- make_seq(CHR2_ord, "force")