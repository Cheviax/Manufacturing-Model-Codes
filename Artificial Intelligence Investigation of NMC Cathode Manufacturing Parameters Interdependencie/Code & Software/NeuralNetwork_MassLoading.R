# Predire avec l'apprentissage automatique

library("e1071")
library(caret)
library(plotly)

setwd("C:\Users\Teo Lombardo\Desktop\Paper_to_submit\ML\ML analysis")
# Charger les données

mass <- read.csv("dataset_mass_loading_low_medium_high.csv")

colnames(mass) <- c("mLoad", "Composition", "StoL", "Viscosity")

# Definir une graine pour rendre le hasard reproductible

set.seed(42)

# Echantillonner au hasard 46 des 58 index de lignes

ind <- sample(
  x = 1:66,
  size=22)

mass44 <- mass[ind,]
mass <- mass44

set.seed(20)


data=mass
data$class=data$mLoad

standardiser <- function(x){
  (x-min(x))/(max(x)-min(x))
}

data[, 2:4] <- lapply(data[, 2:4], standardiser)

data$Low <- c(data$mLoad == "Low")
data$Medium <- c(data$mLoad == "Medium")
data$High <- c(data$mLoad == "High")

formula=as.formula("Low + Medium + High ~ Composition + StoL + Viscosity")

library(caret)
ind=createDataPartition(data$class, times = 1,p = 0.8,list=FALSE)


a=data[ind,]
t=data[-ind,]
dim(a);dim(t)


library(neuralnet)
a$class=NULL
dim(a)
inet <- neuralnet(formula, a,linear.output = TRUE, hidden= c(5,3),lifesign="full",stepmax=10000000)
plot(inet, rep="best")


predict <- compute(inet, t[2:4])
plot(predict$net.result)
dim(t)
n=nrow(t)
result <- which.max(predict$net.result)
for (i in 1:n) { result[i] <- which.max(predict$net.result[i,]) }
for (i in 1:n) { if (result[i]==1) {result[i] = "Low"} }
for (i in 1:n) { if (result[i]==2) {result[i] = "Medium"} }
for (i in 1:n) { if (result[i]==3) {result[i] = "High"} }

table(result,t[,1])

confusionMatrix(as.factor(result),t[,1] )