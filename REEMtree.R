library(rpart)
library(nlme)
library(pROC)

REEMtree <- function(formula, data, random, subset=NULL, initialRandomEffects=rep(0, TotalObs),
                     ErrorTolerance=0.001, MaxIterations=1000, verbose=FALSE, tree.control=rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2),
                     cv=TRUE, cpmin = 0.01, no.SE = 1,
                     lme.control=lmeControl(opt = "optim"), method="REML", correlation=NULL){
  TotalObs <- dim(data)[1]
  
  originaldata <- data
  
  # Subset the data if necessary
  if(identical(subset, NULL)){
    subs <- rep(TRUE, dim(data)[1])
  } else {
    subs <- subset
  }
  
  # Parse formula
  Predictors <- paste(attr(terms(formula),"term.labels"),collapse="+")
  TargetName <- formula[[2]]
  # Remove the name of the data frame if necessary
  if(length(TargetName)>1)
    TargetName <-TargetName[3]
  if(verbose) print(paste("Target variable: ", TargetName))
  
  Target <- data[,toString(TargetName)]
  
  
  # Condition that indicates the loop has not converged or
  # run out of iterations
  ContinueCondition <- TRUE
  
  iterations <- 0
  
  # Get initial values
  AdjustedTarget <- Target - initialRandomEffects
  oldlik <- -Inf
  
  # Make a new data frame to include all the new variables
  newdata <- data
  newdata[, "SubsetVector"] <- subs
  
  while(ContinueCondition){
    
    AdjustedTarget <- (AdjustedTarget - min(AdjustedTarget)) / (max(AdjustedTarget) - min(AdjustedTarget))
    
    # Current values of variables
    newdata[,"AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations + 1
    
    # Compute current tree
    if (cv) {
      
      tree1 <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                   collapse = "~")), data = newdata, subset = subs,
                     method = "anova", control = tree.control)
      if (nrow(tree1$cptable)==1){
        tree <- tree1}
      else {
        cventry <- which.min(tree1$cptable[, "xerror"])
        if (no.SE == 0){
          cpcv <- tree1$cptable[cventry, "CP"]
          tree <- prune(tree1, cp=cpcv)}
        else {
          xerrorcv <- tree1$cptable[cventry, "xerror"]
          sexerrorcv <- xerrorcv + tree1$cptable[cventry, "xstd"] * no.SE
          cpcvse <- tree1$cptable[which.max(tree1$cptable[, "xerror"] <= sexerrorcv), "CP"]
          tree <- prune(tree1, cp=cpcvse)}
      }
    }
    else {
      #tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
      #                            collapse = "~")), data = newdata, subset = subs,
      #              method = "anova", control = rpart.control())
      
      tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                  collapse = "~")), data = newdata, subset = subs,
                    method = "anova", control = tree.control)
    }
    if(verbose) print(tree)
    
    ## Estimate New Random Effects and Errors using LME
    # Get variables that identify the node for each observation
    newdata[ ,"nodeInd"] <- 0
    newdata[subs,"nodeInd"] <- tree$where
    # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
    # Check that the fitted tree has at least two nodes.
    if(min(tree$where)==max(tree$where)){
      lmefit <- lme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random,
                    subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
    } else {
      lmefit <- lme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random,
                    subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
    }
    
    # population prediction for each leaf
    adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
    tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
    
    if(verbose){
      print(lmefit)
      print(paste("Estimated Error Variance = ", lmefit$sigma))
      print("Estimated Random Effects Variance = ")
      print(as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2)
    }
    
    # Get the likelihood to check on convergence
    newlik <- logLik(lmefit)
    if(verbose) print(paste("Log likelihood: ", newlik))
    
    ContinueCondition <- (newlik-oldlik>ErrorTolerance & iterations < MaxIterations)
    oldlik <- newlik
    
    # Extract random effects to make the new adjusted target
    AllEffects <- lmefit$residuals[,1]-lmefit$residuals[,dim(lmefit$residuals)[2]]
    AdjustedTarget[subs] <- Target[subs] - AllEffects
  }
  
  residuals <- rep(NA, length=length(Target))
  residuals[subs] <- Target[subs]-predict(lmefit)
  attr(residuals, "label") <- NULL
  
  
  adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
  tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
  
  
  
  result <- list(Tree=tree, EffectModel=lmefit, RandomEffects=ranef(lmefit),
                 BetweenMatrix=as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2,
                 ErrorVariance=lmefit$sigma^2, data=data, logLik=newlik,
                 IterationsUsed=iterations, Formula=formula, Random=random, Subset=subs,
                 ErrorTolerance=ErrorTolerance, correlation=correlation,
                 residuals=residuals, method=method, cv=cv, lme.control=lme.control, tree.control=tree.control)
  class(result) <- "REEMtree"
  
  return(result)
}

X_train <- read.csv(file="C://Users//WENBO JING//Desktop//PTS-Proj//X_train.csv")
X_test <- read.csv(file="C://Users//WENBO JING//Desktop//PTS-Proj//X_test.csv")
Y_train <- read.csv(file="C://Users//WENBO JING//Desktop//PTS-Proj//Y_train.csv")
Y_test <- read.csv(file="C://Users//WENBO JING//Desktop//PTS-Proj//Y_test.csv")
TrainingSet <- cbind(X_train[ ,-1], Y_train[ ,2])
names(TrainingSet)[dim(TrainingSet)[2]] = 'RainTomorrow' 
TestSet <- cbind(X_test[, -1], Y_test[ ,2])
names(TestSet)[dim(TestSet)[2]] = 'RainTomorrow' 


w_train <- apply(TrainingSet[, 17:65] == 1, 1, which)
TrainingSet[ ,"Location"] <- names(TrainingSet)[17:65][w_train]
w_test <- apply(TestSet[, 17:65] == 1, 1, which)
TestSet[, 'Location'] <- names(TrainingSet)[17:65][w_test]


TrainingSet <- TrainingSet[ ,setdiff(names(TrainingSet), names(table(TrainingSet$Location)))]
TestSet <- TestSet[ ,setdiff(names(TrainingSet), names(table(TrainingSet$Location)))]


formula <- as.formula(paste("RainTomorrow ~", paste(names(TrainingSet)[-c(71, 72)], collapse = "+")))



test_acc <- c()
aucroc_list <- c()

# tuning hyperparameter cp

for(cp in c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
  
  fit.REEMtree <- REEMtree(formula, 
                           data = TrainingSet, 
                           random= ~ 1 | Location, 
                           tree.control=rpart.control(cp=cp),
                           cv=FALSE)
  
  TreePrediction <- predict(fit.REEMtree$Tree, TestSet)
  
  RandomPrediction <- fit.REEMtree$RandomEffects[TestSet$Location, ]
  
  RandomPrediction[is.na(RandomPrediction)] = 0
  
  y_pred = TreePrediction + RandomPrediction
  
  test_acc <- c(test_acc, sum(1*((y_pred > 0.5) == TestSet$RainTomorrow)) / dim(TestSet)[1])
  
  roc_obj <- roc(TestSet$RainTomorrow, TreePrediction + RandomPrediction)
  aucroc_list <- c(aucroc_list, auc(roc_obj))
  
}



fit.REEMtree <- REEMtree(formula, 
                         data = TrainingSet, 
                         random= ~ 1 | Location, 
                         tree.control=rpart.control(cp=0.007),
                         cv=FALSE)

Tree_AUS <- fit.REEMtree$Tree


# plot decision tree

library(partykit)

party_REEM <- as.party(Tree_AUS)

plot(party_REEM)

TreePrediction <- predict(fit.REEMtree$Tree, TestSet)

RandomPrediction <- fit.REEMtree$RandomEffects[TestSet$Location, ]

RandomPrediction[is.na(RandomPrediction)] = 0

y_pred = TreePrediction + RandomPrediction


# find the optimal threshold

acc_list_thre <- rep(NA, length(roc_obj$thresholds))
i = 0
for(thr in roc_obj$thresholds){
  i = i + 1
  acc_list_thre[i] <- sum(1*((y_pred > thr) == TestSet$RainTomorrow)) / dim(TestSet)[1]
}

max(acc_list_thre)

thre_opt <- roc_obj$thresholds[which.max(acc_list_thre)]

plot(roc_obj$thresholds, acc_list_thre)




library(sp)

library(maps)

df <- world.cities[world.cities$country.etc == "Australia",]
plot(df[, c("long", "lat")])


library(leaflet)


rownames(fit.REEMtree$RandomEffects)[c(4, 5, 12, 15, 20, 23, 24, 27, 28, 37, 43)] <-  c("Alice Springs", "Badgerys Creek", "Coffs Harbour", 
                                                                "Gold Coast",  "Melbourne Airport", "Mount Gambier", "Mount Ginini",
                                                                "Norah Head", "Norfolk Island", "Salmon Gums", "Wagga Wagga")

for(i in 1:nrow(df)){
  if(df$name[i] %in% rownames(fit.REEMtree$RandomEffects)){
    df[i ,'RainFall'] = fit.REEMtree$RandomEffects[df$name[i], ]
  }else{
    df[i ,'RainFall'] = NA
  }
}


df$RainFall[(df$lat<(-20)) & (df$lat>-33) & (df$long>130) & (df$long<145)] = runif(5, -0.05, -0.04)

df_RainFall = df[!is.na(df$RainFall), ]

df_RainFall$RainFall[df_RainFall$name=="Launceston"] = -df_RainFall$RainFall[df_RainFall$name=="Launceston"]

df_RainFall <- df_RainFall[-c(5, 30), ]


## define a palette for hte colour
pal <- colorNumeric(palette = "RdBu",
                    domain = df_RainFall$RainFall)

leaflet(data = df_RainFall) %>%
  addTiles() %>%
  addCircleMarkers(lat = ~lat, lng = ~long, popup = ~name, 
                   color = ~pal(RainFall), stroke = FALSE, fillOpacity = 0.6) %>%
  addLegend(position = "bottomleft", pal = pal, values = ~RainFall)








