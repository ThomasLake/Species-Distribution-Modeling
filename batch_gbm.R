#!/usr/bin/env Rscript

#Library loading
cat("Loading Libraries\n")
library(sp)
library(raster)
library(rgdal) 
library(dismo)
library(rJava)
library(rgeos)
library(viridis)
library(maptools)
library(maps)
library(mapdata)
library(httr)
library(ggplot2)
library(prism)
library(gbm)
library(vegan)
library(alphahull)
library(ggplot2)
library(rmaxent)
library(ROCR)
library(vcd)

cat("Loading Environmental Data\n")
#Load in the environmental data
BioMod_all <- brick("c:/users/tom lake/desktop/invasive species/modeling/env_stacks/biomond_1975.grd")

mean.temp.bio1 <- raster(BioMod_all, layer = 1) #(Celcius)
mean.temp.day.range.bio2 <- raster(BioMod_all, layer = 2) #Mean diurnal temperature range (mean(period max-min)) (Celcius)
isotherm.bio3 <- raster(BioMod_all, layer = 3) #Isothermality (Bio2/Bio7)
temp.seas.bio4 <- raster(BioMod_all, layer = 4) #Temperature seasonality (C of V)
max.temp.warm.w.bio5 <- raster(BioMod_all, layer = 5) #Maximum temperature of the warmest week (Celcius)
min.temp.cold.w.bio6 <- raster(BioMod_all, layer = 6) #Minimum temperature of the coldest week (Celcius)
ann.temp.range.bio7 <- raster(BioMod_all, layer = 7) #Annual temperature range (Bio5-Bio6) (Celcius)
mean.temp.wet.q.bio8 <-raster(BioMod_all, layer = 8) #Mean temperature of wettest quarter (Celcius)
mean.temp.dry.q.bio9 <- raster(BioMod_all, layer = 9) #Mean temperature of driest quarter (Celcius)
mean.temp.warm.q.bio10 <- raster(BioMod_all, layer = 10) #Mean temperature of warmest quarter (Celcius)
mean.temp.cold.q.bio11 <- raster(BioMod_all, layer = 11) #Mean temperature of the coldest quarter (Celcius)
ann.precip.bio12 <- raster(BioMod_all, layer = 12) #Annual precipitation (mm)
precip.wet.w.bio13 <- raster(BioMod_all, layer = 13) #Precipitation of the wettest week (mm)
precip.dry.w.bio14 <-raster(BioMod_all, layer = 14) #Preciptitation of the driest week (mm)
precip.seas.bio15 <- raster(BioMod_all, layer = 15) #Preciptitation seasonality
precip.wet.q.bio16 <- raster(BioMod_all, layer = 16) #Precipitation of the wettest quarter (mm)
precip.dry.q.bio17 <- raster(BioMod_all, layer = 17) #Preciptitation of the driest quarter (mm)
precip.warm.q.bio18 <- raster(BioMod_all, layer = 18) #Precipitation of the warmest quarter (mm)
precip.cold.q.bio19 <- raster(BioMod_all, layer = 19) #Precipitation of the coldest quarter (mm)
ann.mean.rad.bio20 <- raster(BioMod_all, layer = 20) #(W m^-2)
hi.w.rad.bio21 <- raster(BioMod_all, layer = 21) #(W m^-2)
low.w.rad.bio22 <- raster(BioMod_all, layer = 22) #(W m^-2)
rad.seas.bio23 <- raster(BioMod_all, layer = 23) #(C of V)
rad.wet.q.bio24 <- raster(BioMod_all, layer = 24) #(W m^-2)
rad.dry.q.bio25 <- raster(BioMod_all, layer = 25) #(W m^-2)
rad.warm.q.bio26 <- raster(BioMod_all, layer = 26) #(W m^-2)
rad.cold.q.bio27 <- raster(BioMod_all, layer = 27) #(W m^-2)
ann.mean.moist.bio28 <- raster(BioMod_all, layer = 28) #Annual mean moisture index
hi.w.moist.bio29 <- raster(BioMod_all, layer = 29) #Highest weekly moisture index
low.w.moist.bio30 <- raster(BioMod_all, layer = 30) #Lowest weekly moisture index
moist.seas.bio31 <- raster(BioMod_all, layer = 31) #Moisture index seasonality (C of V)
mean.moist.wet.q.bio32 <- raster(BioMod_all, layer = 32) #Mean moisture index of the wettest quarter
mean.moist.dry.q.bio33 <- raster(BioMod_all, layer = 33) #Mean moisture index of the driest quarter
mean.moist.warm.q.bio34 <- raster(BioMod_all, layer = 34) #Mean moisture index of the warmest quarter
mean.moist.cold.q.bio35 <- raster(BioMod_all, layer = 35) #Mean moisture index of the coldest quarter
BioMod.PC1.bio36 <- raster(BioMod_all, layer = 36) #Principle component 1
BioMod.PC2.bio37 <- raster(BioMod_all, layer = 37) #Principle component 2
BioMod.PC3.bio38 <- raster(BioMod_all, layer = 38) #Principle component 3
BioMod.PC4.bio39 <- raster(BioMod_all, layer = 39) #Principle component 4
BioMod.PC5.bio40 <- raster(BioMod_all, layer = 40) #Principle component 5

BioMod_gbm <- stack(mean.temp.bio1, mean.temp.day.range.bio2, isotherm.bio3, temp.seas.bio4, max.temp.warm.w.bio5, min.temp.cold.w.bio6, ann.temp.range.bio7, mean.temp.wet.q.bio8, mean.temp.dry.q.bio9, mean.temp.warm.q.bio10, mean.temp.cold.q.bio11, ann.precip.bio12, precip.seas.bio15, precip.wet.q.bio16, precip.dry.q.bio17, precip.warm.q.bio18, precip.cold.q.bio19, ann.mean.rad.bio20, rad.seas.bio23, rad.wet.q.bio24, rad.dry.q.bio25, rad.warm.q.bio26, rad.cold.q.bio27, ann.mean.moist.bio28, moist.seas.bio31, mean.moist.wet.q.bio32, mean.moist.dry.q.bio33, mean.moist.warm.q.bio34, mean.moist.cold.q.bio35)

#Load GBM Species Dataset from generate_gbm_dataset.R output
pa.gbm.dataset <- read.csv("c:/users/tom lake/desktop/maxentoutputs/model1/PA_GBM_FullDataset_wEnvData.csv")
coordinates(pa.gbm.dataset) <- ~ lon + lat
projection(pa.gbm.dataset) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#function to perform gbm simplify
gbmSimplify <- function(gbmstep){
  
  simpl <- gbm.simplify(gbmstep)
  
  return(simpl)
}

#core function to perform gbm modeling
gbmStep <- function(dataset, env, arguments){
  
  len <- NCOL(dataset) - 1
  tc <- arguments[[1]]
  lr <- arguments[[2]]
  ntrees <- arguments[[3]]
  bf <- arguments[[4]]
  
  step_1 <- gbm.step(data = dataset, gbm.x = 3:len, gbm.y = NCOL(dataset), family = "bernoulli", tree.complexity = tc, learning.rate = lr, bag.fraction = bf)
  simplified <- gbmSimplify(step_1)
  step_2 <- gbm.step(data = dataset, gbm.x = simplified$pred.list[[11]], gbm.y = NCOL(dataset), tree.complexity = tc, learning.rate = lr, bag.fraction = bf)
  
  #predictions
  pred_1.1 <- dismo::predict(env, step_1, n.trees = ntrees, type = "response")
  pred_1.2 <- dismo::predict(env, step_1, n.trees = step_1$gbm.call$best.trees, type = "response")
  
  pred_2.1 <- dismo::predict(env, step_2, n.trees = ntrees, type = "response")
  pred_2.2 <- dismo::predict(env, step_2, n.trees = step_2$gbm.call$best.trees, type = "response")
  
  predictions <- stack(pred_1.1, pred_1.2, pred_2.1, pred_2.2)
  return(list(step_1, step_2, predictions))
}

#function to specify gbm arguments
gbmArgs <- function(treec, learnr, ntree, bagfrac){
  
  #tree complexity
  tc = treec
  #learning rate
  lr = learnr
  #number of trees
  ntrees = ntree
  #bag fraction
  bf = bagfrac

  args = list(tc, lr, ntrees, bf)
  
  return(args)
}

#function calculates model evaluation and raster statistics
calculateStatistics <- function(model, dataset, n){
  
  cat('calculating model statistics')
  eval_stats <- data.frame(matrix(ncol=4))
  colnames(eval_stats) <- c("auc", "kappa", "tss", "ntrees", "learning rate")
  
  for(j in 1:length(model)){
    
    #Calculate auc, tss, and kappa statistics
    step1model <- model[[j]][[n]]
    ntree <- step1model$gbm.call$best.trees
    print(ntree)
    preds <- predict(step1model, dataset, n.trees = step1model$gbm.call$best.trees, type="response")
    c <- calc.deviance(obs=dataset$presence, pred = preds, calc.mean=TRUE)
    d <- cbind(dataset$presence, preds)
    trainpp <- d[d[,1]==1, 2]
    print(trainpp)
    bb <- d[d[,1]==0, 2]

    #auc
    e <- evaluate(p=trainpp, a=bb)
    auc <- e@auc
    
    #evaluation methods
    combined <- c(trainpp, bb)
    label <- c(rep(1,length(trainpp)), rep(0,length(bb)))
    pred <- prediction(combined,label)
    perf <- performance(pred, "tpr", "fpr")
    fpr = perf@x.values[[1]]
    tpr = perf@y.values[[1]]
    sum = tpr + (1-fpr)
    index = which.max(sum)
    thresh = perf@alpha.values[[1]][[index]]
    
    #confusion matrix
    con <- cbind(c(length(trainpp[trainpp>=thresh]), length(trainpp[trainpp<thresh])),
                       c(length(bb[bb>=thresh]), length(bb[bb<thresh])))
    #kappa
    kappavalue <- Kappa(con)
    kappa <- kappavalue$Unweighted[[1]]
    
    #tss
    sensitivity<-con[1,1]/(con[1,1]+con[2,1])
    specificity<-con[2,2]/(con[2,2]+con[2,1])
    tss<-sensitivity+specificity-1
    
    eval_stats[j,] <- c(auc, kappa, tss, ntree)
    
  }
  return(eval_stats)
}

#function for raster statistics
gbmStatistics <- function(gbmModels, dataset){
  
  #raster statistics
  cat('calculating raster statistics')
  predictions <- stack()
  gbmMean <- stack()
  gbmRange <- stack()
  
  for(i in 1:length(gbmModels)){
    predictions <- stack(gbmModels[[i]][[3]])
  }
  
  gbmMean <- calc(predictions, mean)
  plot(gbmMean)
  gbmRange <- calc(predictions, range)
  plot(gbmRange)
  
  #number of trees
  
  #model evaulation statistics
  eval_stats <- calculateStatistics(gbmModels, dataset, 1)
  
  
return(eval_stats)
}

#function for batch gbm models
gbmBatch <- function(dataset, env, n, tc, lr, ntrees, bf){
  
  #default arguments
  arguments <- gbmArgs(tc, lr, ntrees, bf)
  #default model list
  models <- list()

  for(i in 1:n){
  
    gbmModel <- gbmStep(dataset, env, arguments)
  
    models[[i]] <- gbmModel
    
    lr = lr + 0.001
    ntrees = ntrees + 50
    arguments <- gbmArgs(tc, lr, ntrees, bf)
    
    if(i == n){
      stats <- gbmStatistics(models, dataset)
      print(stats)
    }
    
  }
  
  return(models)
}


###Ideas
#gbmstep is core function that performs operations [x]
#need to add functionality around gbmstep [x]
#input varying datasets [future]
#ndrops & number of variables [x]
#input tc and lr and bf [x]
#output best number of trees created for each learning rate [x]
#relationship between learning rate and n trees
#number of trees in model statistics output
#integrate plotting/drawing function into maxent batch and maybe gbm batch
#look into msi biomod2 debugging
#methods for batch scripts
#general methods for rob
#R methods for amy


###Testing

batchtest <- gbmBatch(pa.gbm.dataset, BioMod_gbm, 2, tc = 5, lr = 0.005, ntrees = 5000, bf = .5)

stats <- gbmStatistics(batchtest, pa.gbm.dataset)


batchtest[[1]]
#batchtest[[1]][[1]] is gbm first step call
#batchtest[[1]][[2]] is gbm second step call
#batchtest[[1]][[3]] is gbm prediction raster stack


b <- batchtest[[1]][[1]]

preds <- predict(pa.tc5.l004.simpl.1,pa.gbm.dataset, n.trees = pa.tc5.l004.simpl.1$n.trees, type="response")
c <- calc.deviance(obs=pa.gbm.dataset$presence, pred = preds, calc.mean=TRUE)
d <- cbind(pa.gbm.dataset$presence, preds)
trainpp <- d[d[,1]==1, 2]
bb <- d[d[,1]==0, 2]

#auc
cat('calculating aic')
e <- evaluate(p=trainpp, a=bb)
auc <- e@auc

#evaluation methods
combined <- c(trainpp, bb)
label <- c(rep(1,length(trainpp)), rep(0,length(bb)))
pred <- prediction(combined,label)
perf <- performance(pred, "tpr", "fpr")
fpr = perf@x.values[[1]]
tpr = perf@y.values[[1]]
sum = tpr + (1-fpr)
index = which.max(sum)
thresh = perf@alpha.values[[1]][[index]]

#confusion matrix
confusion <- cbind(c(length(trainpp[trainpp>=thresh]), length(trainpp[trainpp<thresh])),
                   c(length(bb[bb>=thresh]), length(bb[bb<thresh])))
#kappa
cat('calculating kappa')
kappavalue <- Kappa(confusion)
kappa <- kappavalue$Unweighted[[1]]

#tss
cat('calculating tss')
con <- confusion
sensitivity<-con[1,1]/(con[1,1]+con[2,1])
specificity<-con[2,2]/(con[2,2]+con[2,1])
tss<-sensitivity+specificity-1


eval_stats[j,] <- c(auc, kappa, tss)


#plot predictions test
ics <- ic(p2$prediction_raw, train_p, me.pa.BioModPA5.100.std.tr1.b1)

Stats.eval.names <- c("AUC-test", "Threshold", "Kappa", "TSS", "AICc", "BIC")
Stats.eval.values <- rbind(round(aucvalue, digits = 2), round(cutoff, digits =2), round(kappavalue$Unweighted[1], digits = 2), round(tssvalue, digits = 2), round(ics[5], digits = 2), round(ics[6], digits = 2))

Stats.eval.dataframe <- data.frame(matrix(nrow = length(Stats.eval.names), ncol = 2))
Stats.eval.dataframe[1] <- Stats.eval.names
Stats.eval.dataframe[2] <- Stats.eval.values
colnames(Stats.eval.dataframe) <- c("Statistic", "Value")
stat.eval.table <- ggtexttable(Stats.eval.dataframe, rows = NULL)


#Plot maxent model
#levelplot(p2$prediction_cloglog, margin = FALSE, maxpixels = 1e5, col.regions = function(x) rev(terrain.colors(x)), at = seq(0,1, len = 256), main = "Palmer's Amaranth\nModel Built with Current Distribution\n Train = 1; Beta = 1", sub = textGrob(paste("AUC = ", round(aucvalue, digits = 2), "\nKAPPA = ", round(kappavalue$Unweighted[1], digits = 2), "\nTSS = ", round(tssvalue, digits = 2), "\nAICc = ", round(ics[5], digits = 2), "\nBIC = ", round(ics[6], digits = 2)), x = unit(0, "npc"), y = unit(0.5, "npc"), just = "left", rot = 0, gp = gpar(cex = 0.5))) + spplot(states.lines.crop, col.regions = "gray20") + spplot(subsamp.pa.all.coords.biomod, pch = 19, cex = 0.1, col.regions = "red")

#Modifying the key to put in model evals and parameters
#evals <- list(title = "Model Evals", space = "bottom", columns = 2, text = list(paste("AUC = ", round(aucvalue, digits = 2), "\nKAPPA = ", round(kappavalue$Unweighted[1], digits = 2), "\nTSS = ", round(tssvalue, digits = 2), "\nAICc = ", round(ics[5], digits = 2), "\nBIC = ", round(ics[6], digits = 2))), cex.title = 1, cex = 0.7)

base.prediction.plot <-levelplot(p2$prediction_cloglog, margin = FALSE, maxpixels = 1e5, col.regions = function(x) rev(terrain.colors(x)), at = seq(0,1, len = 256), main = "Palmer's Amaranth")
state.line.plot <- spplot(states.lines.crop, col.regions = "gray20")
occurrence.pt.plot <- spplot(subsamp.pa.all.coords.biomod, pch = 19, cex = 0.1, col.regions = "red")
prediction.plot <- base.prediction.plot + state.line.plot + occurrence.pt.plot
text <- paste("Model Built with Current Distribution", "Train = 1; Beta = 1", sep = " ")
text.p <- ggparagraph(text = text, face = "italic", color = "black")

cols <- 1
rows <- 4

figure <- multi_panel_figure(width = 135, columns = cols, height = 100, rows = rows)
(figure %<>% fill_panel())
(figure %<>% fill_panel(text.p))



#ggarrange(prediction.plot, text.p, n.col = 1, nrow = 2, heights = c(1, 0.5)) #Does not work because the map is a trellis object


ztable(Clim.import.dataframe)
ztable(Stat.eval.dataframe)







