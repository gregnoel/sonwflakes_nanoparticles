setwd("C://Users/Abeille/Downloads/data2")

library(tidyverse)
library(caret)
library(MASS)
library(ggplot2)






###Other test
library(klaR)
library(psych)
library(MASS)
library(ggord)
library(devtools)


#Import and dataset for peptide selection with all the peptide
Five_prot<-read.csv2("Matrix_Training_5Target.csv")
Five_prot$Protein<-as.factor(Five_prot$Protein)
str(Five_prot)
Five_prot<-Five_prot[-1]

#Graphique ggord
FiveProtLDA <- lda(Protein~., Five_prot)
FiveProtLDA
kk<-ggord(FiveProtLDA, Five_prot$Protein,  ellipse_pro = 0.95, arrow = NULL, txt = NULL)
kk
# ylim = c(-12, 10),xlim = c(-27,15), arrow = NULL, txt = NULL



#With only the Covid, Flue and Healthy
Moko<-read.csv2("Matrix_Training_3Target.csv")
Moko$Protein<-as.factor(Moko$Protein)
str(Moko)
Moko<-Moko[-1]
#Moko<-Moko[-c(7),]

#Graphic from LDA for the 3 diagnosis
linear <- lda(Protein~., Moko)
linear
kk<-ggord(linear, Moko$Protein, ellipse_pro = 0.95, arrow = NULL, txt = NULL)
kk

#ggsave("Graph_LDA_Trained.png", kk, width = 12, height = 8, dpi = 750 , units = "cm")


######Preparation of tested matrix
MatrixTest<-read.csv2("MatrixTest.csv")
MatrixTest$Protein<-as.factor(MatrixTest$Protein)
str(MatrixTest)
MatrixTest<-MatrixTest[-1]
######
#Assignation of training matrix and matrix to test

train.data <- Moko
test.data <- MatrixTest


# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)




# Fit the model
model <- lda(Protein~., data = train.transformed)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$Protein)


model <- lda(Protein~., data = train.transformed)
model


plot(model)

predictions <- model %>% predict(test.transformed)
names(predictions)



# Predicted classes
head(predictions$class, 9)
# Predicted probabilities of class memebership.
head(predictions$posterior, 9) 
# Linear discriminants
head(predictions$x, 9) 

Prediction<-as.data.frame(predict(model)$x)
lda.data <- cbind(train.transformed, Prediction)
lda.data2<-cbind(test.transformed, predictions$x)
#write.csv2(lda.data,"lda_data_train.csv")
#write.csv2(lda.data2, "lda_data_test.csv")

total_lda<-read.csv2("lda_data_total.csv")
total_lda$Protein<-as.factor(total_lda$Protein)


kp<-ggplot(total_lda, aes(LD1, LD2,color = Protein)) +
  geom_point(size = 3) +  stat_ellipse(level = 0.99)  + labs(x = "LD1 (87.45%)", y = "LD2 (12.55%)")

kp<- kp + theme_bw(base_size = 12)
kp



#####With blood proteins
EightProt<-read.csv2("Matrix_8Prot_Control_target.csv")
EightProt$Protein<-as.factor(EightProt$Protein)
str(EightProt)


#Graphic from LDA for the 3 diagnosis
linear <- lda(Protein~., EightProt)
linear
kk<-ggord(linear, EightProt$Protein, ellipse_pro = 0.75, arrow = NULL, txt = NULL)
kk

#ggsave("Graph_LDA_Trained.png", kk, width = 12, height = 8, dpi = 750 , units = "cm")




###Preparation from Borcard 2018 LDA


# Load the required packages
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)

spe.norm<-decostand(Moko[-1],"rank")

gr <- Moko$Protein

prote<-as.matrix(Moko[-1])
# Verify multivariate homogeneity of within-group covariance
# matrices using the betadisper() function {vegan}
prote.d1 <- dist(spe.norm, method = "manhattan")
(env.MHV <- betadisper(prote.d1, Moko$Protein))
permutest(env.MHV)

Wilks.test(prote, gr)


# Computation of LDA - identification functions (on unstandardized
# variables)
env.pars3.df <- as.data.frame(prote)
(spe.lda <- lda(Moko$Protein ~ ., data = env.pars3.df))
# The result object contains the information necessary to interpret
# the LDA
summary(spe.lda)
# Display the group means for the 3 variables
spe.lda$means
# Extract the unstandardized identification functions (matrix C,
# eq. 11.33 in Legendre and Legendre 2012)
(C <- spe.lda$scaling)
# Classification of two new objects (identification)#prediction

(predict.new <- predict(spe.lda, newdata = MatrixTest))


# Computation of LDA - discrimination functions (on standardized
# variables)
env.pars3.sc <- as.data.frame(scale(env.pars3.df))
spe.lda2 <- lda(gr ~ ., data = env.pars3.sc)
# Display the group means for the 3 variables
spe.lda2$means
# Extract the classification functions
(C2 <- spe.lda2$scaling)
# Compute the canonical eigenvalues
spe.lda2$svd^2
# Position the objects in the space of the canonical variates
(Fp2 <- predict(spe.lda2)$x)
# Classification of the objects
(spe.class2 <- predict(spe.lda2)$class)
# Posterior probabilities of the objects to belong to the groups
# (rounded for easier interpretation)
(spe.post2 <- round(predict(spe.lda2)$posterior, 2))
# Contingency table of prior versus predicted classifications
(spe.table2 <- table(gr, spe.class2))
# Proportion of correct classification (classification success)
diag(prop.table(spe.table2, 1))


# plot.lda
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/plot.lda.R')


# Plot the LDA results using the homemade function plot.lda()
plot.lda(lda.out = spe.lda2,
         groups = gr,
         plot.sites = 2,
         plot.centroids = 1,
         mul.coef = 2.35
)

# LDA with jackknife-based classification (i.e., leave-one-out
# cross-validation)
(spe.lda.jac <-
    lda(gr ~ .,
        data = env.pars3.sc,
        CV = TRUE))
summary(spe.lda.jac)
# Numbers and proportions of correct classification
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
# Classification success
diag(prop.table(spe.jac.table, 1))









library(e1071)

confusionMatrix(Protein~., predict(linear, Moko))












pairs.panels(Moko[2:9],
             gap = 0,
             bg = c("red", "green", "blue","purple","gold","black","brown","orange")[Moko$Protein],
             pch = 21)

set.seed(123)
ind <- sample(2, nrow(Moko),
              replace = TRUE,
              prob = c(1.00, 0))
training <- Moko[ind==1,]
testing <- Moko[ind==2,]



linear <- lda(Protein~., Moko)
linear


p <- predict(linear, training)
ldahist(data = p$x[,1], g = training$Protein)

ldahist(data = p$x[,2], g = training$Protein)
ggord(linear, Moko$Protein, ylim = c(-10, 10),xlim = c(-20,20), ellipse_pro = 0.80, shape = Moko$Protein)

#partimat(Protein~., data = training, method = "lda")

p1 <- predict(linear, training)$class
tab <- table(Predicted = p1, Actual = training$Protein)
tab

sum(diag(tab))/sum(tab)


p2 <- predict(linear, testing)$class
tab1 <- table(Predicted = p2, Actual = testing$Protein)
tab1
sum(diag(tab1))/sum(tab1)




###PCA

install.packages(c("FactoMineR", "factoextra"))


library("FactoMineR")
library("factoextra")

Moko_t<-na.omit(Moko)
res.pca<-PCA(Moko[,-1], graph = T)


eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_pca_var(res.pca, col.var = "black")

fviz_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Moko$Protein, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#00bba5","#00bb13","#bb7000","#bb00a5","#bb0019","#0041bb" ),
             addEllipses = T, # Concentration ellipses
             legend.title = "Protein"
)


###Retest

install_github("Displayr/flipMultivariates")

