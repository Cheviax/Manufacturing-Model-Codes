# Predire avec l'apprentissage automatique

setwd("C:\Users\Teo Lombardo\Desktop\Paper_to_submit\ML\ML analysis")
# Charger les données

mass <- read.csv("dataset_mass_loading_low_medium_high.csv")

colnames(mass) <- c("mLoad", "Composition", "StoL", "Viscosity")
# Definir une graine pour rendre le hasard reproductible

set.seed(93)

# Standardiser le dataset

# standardiser <- function(x){
#   (x-min(x))/(max(x)-min(x))
# }
# 
# mass[, 2:4] <- lapply(mass[, 2:4], standardiser)


# Echantillonner au hasard 46 des 58 index de lignes
indexes <- sample(
  x = 1:66,
  size=53)

# Visualiser les index aleatoires

print(indexes)

# Creer un ensemble de formation à partir d'index

train <- mass[indexes,]

# Creer un ensemble de test à partir des index restants

test <- mass[-indexes,]

# Charger le package de l'arbre de décision


library(tree)

# Former un modele d'arbre de décision

model <- tree(
  formula = mLoad ~ Composition+Viscosity+StoL,
  data = train,
  control = tree.control(nobs=53,mincut = 2, minsize = 26, mindev = 0.01),
  method = "class",
  split = c("deviance"))

# Visualiser le modèle

summary(model)


# Visualiser le modele de l'abre de decision

plot(model)


text(model)



# Sauvegarder le modèle de l'arbre

#save(model,file = "Tree.Rdata")

# Sauvegarder les données d'entrainement

#save(train, file = "Train.Rdata")

# Sauvegarder les données complètes

#save(mass, file = "Mass.Rdata")

# Charger la bibliotheque de brasseurs de couleurs
library(RColorBrewer)

# Créer une palette de couleurs
palette <- brewer.pal(3,"Set1")

# Créer un nuage de points coloré par espace

plot(
  x = subset(train$Viscosity, train$Composition==96),
  y = subset(train$StoL, train$Composition==96),
  pch = 19, 
  col = palette[as.numeric(train$mLoad) ],
  main = "",
  xlab = "Viscosity (Pa.Hz)",
  ylab = "S-to-L (%)")

# Tracer les limites de décision

partition.tree(
  tree=model,
  label="mLoad",
  add=TRUE)


# Prédire avec le modèle

predictions <- predict(
  object = model,
  newdata = test,
  type = "class")

# Créer une matrice de confusion

table(
  x = predictions, 
  y = test$mLoad)

# Installer le package caret


library(caret)

# Evaluer les résultats de la prédiction


library(e1071)


confusionMatrix(
  data = predictions,
  reference = test$mLoad)