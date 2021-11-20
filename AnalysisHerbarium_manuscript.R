# set contrasts for session (always do this)
#important
options(contrasts = c("contr.sum", "contr.poly")) 
#install.packages("DHARMa")
library(DHARMa)
library(car)
library(MASS)
library(cowplot)
library(jtools)
library(interactions)
R2logit<- function(model){
  R2<- 1-(model$deviance/model$null.deviance)
  return(R2)
}

#read in trait data
data<-read.table ("~/Dropbox/Documents/ragweed/210721 Herbarium Analysis - Northern Hemisphere.txt", sep="\t", header=T)
Europe<-data[data$Continent=="EUROPE",]
NAm<-data[data$Continent=="NORTH_AMERICA",]


full <- glm(Mature_Male_Inflorescence ~  Day* Latitude*Year, data = Europe, family = "binomial")
summary(full)
Anova(full, type = 3, test="F")
full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Latitude+Day:Year+Latitude:Year +Day:Latitude:Year, data = Europe, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Latitude+Day:Year+Latitude:Year , data = Europe, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Year+Latitude:Year , data = Europe, family = "binomial")
Anova(full, type = 3, test="F")
testDispersion(full)
simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)

#best plot
A<-interact_plot(full, pred=Latitude , modx=Year, plot.points = TRUE, interval=TRUE , y.label = "Probability of Flowers")
B<-interact_plot(full, pred=Day , modx=Year, plot.points = TRUE, interval=TRUE , y.label = "Probability of Flowers")
plot_grid(A, B, labels = c('A', 'B'), label_size = 12)
R2logit(full)

full <- glm(Presence_of_Fruit ~  Day*Latitude *Year, data = Europe, family = "binomial")
summary(full)
Anova(full, type = 3, test="F")
testDispersion(full)
simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year + Day:Latitude +Day:Year + Latitude:Year, data = Europe, family = "binomial")
Anova(full, type = 3, test="F")
testDispersion(full)
simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year + Day:Latitude +Day:Year, data = Europe, family = "binomial")
Anova(full, type = 3, test="F")
R2logit(full)
testDispersion(full)
simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)

#best plot
C<-interact_plot(full, pred=Day , modx=Year, plot.points = TRUE,  interval=TRUE , y.label = "Probability of Fruit")
D<-interact_plot(full, pred=Day , modx=Latitude, plot.points = TRUE,  interval=TRUE , y.label = "Probability of Fruit")

plot_grid(A, B, C, D, labels = c('A', 'B', 'C', 'D'), label_size = 12)

full <- glm(Mature_Male_Inflorescence ~  Day* Latitude*Year, data = NAm, family = "binomial")
summary(full)

full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Latitude+Day:Year+Latitude:Year +Day:Latitude:Year, data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Latitude+Day:Year+Latitude:Year , data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Mature_Male_Inflorescence ~  Day+ Latitude+Year+Day:Year+Day:Latitude , data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
R2logit(full)
testDispersion(full)
simulationOutput <- simulateResiduals(fittedModel = full, plot = F)
plot(simulationOutput)

#best plot
A<-interact_plot(full, pred=Day , modx=Year, plot.points = TRUE, interval=TRUE , y.label = "Probability of Flowers")
B<-interact_plot(full, pred=Day , modx=Latitude, plot.points = TRUE, interval=TRUE , y.label = "Probability of Flowers")

plot_grid(A, B, labels = c('A', 'B'), label_size = 12)

#nothing significant
full <- glm(Presence_of_Fruit ~  Day*Latitude *Year, data = NAm, family = "binomial")
summary(full)
Anova(full, type = 3, test="F")
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year + Day:Latitude +Day:Year + Latitude:Year, data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year + Day:Latitude +Day:Year, data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year + Day:Latitude , data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
full <- glm(Presence_of_Fruit ~  Day+Latitude +Year , data = NAm, family = "binomial")
Anova(full, type = 3, test="F")
R2logit(full)

