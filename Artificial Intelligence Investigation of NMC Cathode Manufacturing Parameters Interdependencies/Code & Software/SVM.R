# Installer le package shiny
library(shiny)

# Installer les librairies utiles pour l'analyse de données
library("e1071")
library(caret)
library(plotly)

# Définir le répertoire de travail

setwd("C:/Users/Teo Lombardo/Desktop/Paper_to_submit/ML/ML analysis")

mass <- read.csv("V2_dataset_mass_loading_low_medium_high.csv")

colnames(mass) <- c("mLoad", "Composition", "StoL", "Viscosity")

# Definir une graine pour rendre le hasard reproductible

set.seed(16)

# Echantillonner au hasard 46 des 58 index de lignes

indexes <- sample(
  x = 1:82,
  size=82*0.8)

# Creer un ensemble de formation à partir d'index

train <- mass[indexes,]

# Creer un ensemble de test à partir des index restants

test <- mass[-indexes,]

# Attach the Data
attach(train)

# Divide data to x (containt the all features) and y only the classes
x <- subset(train, select=-mLoad)
y <- mLoad

# Create SVM Model and show summary

svm_model <- svm(mLoad ~ Composition+StoL+Viscosity,kernel="radial", data=train)
summary(svm_model)

# Create SVM Model and show summary

svm_model1 <- svm(x,y)
summary(svm_model1)

# Run Prediction and you can measuring the execution time in R

pred <- predict(svm_model1,x)
system.time(pred <- predict(svm_model1,x))

# See the confusion matrix result of prediction, using command table to compare the result of SVM prediction and the class data in y variable.
table(pred,y)

# Tuning SVM to find the best cost and gamma
svm_tune <- tune(svm, train.x=x, train.y=y, 
                 kernel="radial", ranges=list(cost=10^(-1:2), gamma=seq(0,1,0.1)))

print(svm_tune)
BestCost <- svm_tune$best.parameters[,c(1)]
BestGamma <- svm_tune$best.parameters[,c(2)]

# After you find the best cost and gamma, you can create svm model again and try to run again
svm_model_after_tune <- svm(mLoad ~ Composition+StoL+Viscosity, data=train, kernel="radial", cost=50, gamma=0.1)
summary(svm_model_after_tune)


# Run Prediction again with new model
pred <- predict(svm_model_after_tune,x)
system.time(predict(svm_model_after_tune,x))


# See the confusion matrix result of prediction, using command table to compare the result of SVM prediction and the class data in y variable
table(pred,y)

# Charger la library de couleur

library(RColorBrewer)

# Créer une palette de couleurs 

palette <- brewer.pal(3,"Set2")

# Créer un code d'interface utilisateur

InterfaceUtil <- fluidPage(
  titlePanel("Mass Loading Level Predictor"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        inputId = "composition",
        label = "Composition (A.M.%)",
        min = 92.7,
        max = 96,
        step = 0.1,
        value = 92.7),
      sliderInput(
        inputId = "viscosity",
        label = "Viscosity (Pas)",
        min = 1.36,
        max = 11.49,
        step = 0.01,
        value = 1.36),
      sliderInput(
        inputId = "StoL",
        label = "Solid-to-liquid Ratio (%)",
        min = 55,
        max = 80,
        step = 0.01,
        value = 65)),
    mainPanel(
      textOutput(
        outputId = "text"),
      plotOutput(
        outputId = "plot"))))

# Créer un code serveur

serveur1 <- function(input,output){
  
  output$text = renderText({
    
    # Créer des prédicteurs
    
    predictors <- data.frame(
      Composition = input$composition,
      Viscosity = input$viscosity,
      StoL = input$StoL)
    
    # Faire des prédictions
    
    prediction = predict(
      object = svm_model_after_tune,
      newdata = predictors,
      type = "class")
    
    
    # Créer un texte de prédiction
    
    paste("The predicted Mass Loading Level is >>>>>>",
          as.character(prediction),"<<<<<<")
    
    
  })
  
  # Charger la bibliotheque de brasseurs de couleurs
  library(RColorBrewer)
  
  # Créer une palette de couleurs
  palette <- brewer.pal(3,"Paired")
  
  output$plot = renderPlot({
    
    plot(svm_model_after_tune, data=train, StoL ~ Viscosity, 
         fill=TRUE, grid = 200,
         slice = list(Composition = input$composition),
         svSymbol ="",dataSymbol="")
    
    # Dessine un prédicteur sur la parcelle
    
    # points(
    #   x = input$viscosity,
    #   y = input$StoL,
    #   col = "red",
    #   pch = 4,
    #   cex = 2,
    #   lwd = 2)
    
    
  })
}

# Créer une application shiny

shinyApp(
  ui = InterfaceUtil,
  
  server = serveur1)