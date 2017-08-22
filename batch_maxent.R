#Maxent Batch Modeling Script
#Thomas Lake & Ryan Briscoe Runquist
#8/8/17

#Libraries
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
library(rmaxent)
library(biomod2)
library(ROCR)
library(vcd)
library(rasterVis)
library(lattice)
library(devtools)
library(dismotools)

#Load Environmental Dataset
super_stack <- brick('c:/users/tom lake/desktop/invasive species/modeling/env_stacks/climond_min_pop_veg_prism.grd')

precip <- raster(super_stack, layer=49)
mean.temp <- raster(super_stack, layer=50)
min.temp <- raster(super_stack, layer=51)
max.temp <- raster(super_stack, layer=52)
dew.point <- raster(super_stack, layer=53)
vp.max <- raster(super_stack, layer=54)
vp.min <- raster(super_stack, layer=55)
Calcium <- raster(super_stack, layer=41)
Potassium <- raster(super_stack, layer=42)
Magnesium <- raster(super_stack, layer=43)
Phosphorous <- raster(super_stack, layer=44)
Organic_Carbon <- raster(super_stack, layer=46)
veg <- raster(super_stack, layer=48)

BioMod_reduced <- stack(precip, mean.temp, dew.point, Phosphorous, veg)

#Load Species Coordinates Dataset
subsamp.pa.all.coords.biomod <- read.csv("c:/users/tom lake/desktop/invasive species/modeling/palmers_coords/palmers_expanded_range_prism.csv", header = TRUE, sep = ",")
coordinates(subsamp.pa.all.coords.biomod) <- ~ lon + lat
projection(subsamp.pa.all.coords.biomod) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")


ModelOutputPath <- function(){
  #Function to set directory for all model outputs
  #User needs to change path manually
  opath = "c:/users/tom lake/desktop/maxentoutputs/"
  return(opath)
}

sampleBackground <- function(spp){
  #Function inputs species coordinates
  #Samples 10,000 background points within 'doughnut'circle
  
  #Draw pseudoabsence sampling cirlces
  cat('Drawing pseudoabsence sampling circles\n')
  circ.big <- circles(spp, d=100000, lonlat=TRUE)
  #circ.little <- circles(spp, d=5000, lonlat=TRUE)
  pol.big <- polygons(circ.big)
  #pol.little <- polygons(circ.little)
  #doughnut <- gDifference(pol.big, pol.little) #makes a new polygon that has 100km circles but excludes the 5km around the species occurrence record
  
  cat('Sampling pseudoabsences\n')
  pseudo.pts <- spsample(pol.big, 10000, type='stratified', iter = 1000)
  projection(pseudo.pts) <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  
  return(pseudo.pts)
}

extractClimateData <- function(spp, env, pseudo){
  #Function extracts climate data for presence & pseudoabsence points
  
  cat('Extracting climate data\n')
  pa.biomod <- extract(env, spp)
  pa.biomod.pseudo <- extract(env, pseudo)
  
  pa.biomod <- data.frame(spp, pa.biomod)

  cat('Removing presence NA from dataframes\n')
  pa.biomod <- na.omit(pa.biomod)
 
  cat('Removing pseudoabsence NA from dataframes\n')
  pa.biomod.pseudo.pts <- pseudo@coords
  colnames(pa.biomod.pseudo.pts) <- c("lon", "lat")
  pa.biomod.pseudo <- data.frame(cbind(pa.biomod.pseudo.pts, pa.biomod.pseudo))
  pa.biomod.pseudo <- na.omit(pa.biomod.pseudo)
  
  pa <- list(pa.biomod, pa.biomod.pseudo)
  return(pa)
}

crossValidation <- function(presence, pseudo, n){
  #Function performs n-fold cross validation, returns lon, lat cross validation list
  cat(sprintf('Performing %s fold cross-validation\n', n))
  group_p <- kfold(presence, n)
  group_a <- kfold(pseudo, n)
  
  CV <- list()
  
  for(i in 1:n){
    
    cv <- i 
    train_p <- presence[group_p!=i, c("lon", "lat")]
    train_a <- pseudo[group_a!=i,c("lon", "lat")]
    test_p <- presence[group_p==i,c("lon", "lat")]
    test_a <- pseudo[group_a==i,c("lon", "lat")]
  
    CV[[i]] <- c(train_p, train_a, test_p, test_a) 
  }
  return(CV) #CV[[1]][1:2] returns [[cross validation 1]][train_p(lon, lat)]
}

maxentArgs <- function(){
  #function to set maxent model arguments
  args <- c('betamultiplier=1','writebackgroundpredictions','responsecurves',
    'jackknife','randomtestpoints=20','replicates=5','replicatetype=bootstrap',
    'outputgrids=FALSE','randomseed=TRUE', 'threads=4')
  return(c(args))
}

runMaxent <- function(presence, pseudo, env, opath, arguments){
  #Function to run dismo MaxEnt models
  maxentModel <- dismo::maxent(env, p = presence, a = pseudo, 
                        path = opath, args=arguments)
  return(maxentModel)
}

maxentStatistics <- function(modelList, env, presences){
  #Function to calculate maxent model statistics
  #Returns data frame of auc, kappa, and tss statistics
  
  cat('Calculating model statistics\n')
  folders <- list.files(ModelOutputPath())
  
  eval_stats <- data.frame(matrix(ncol=3))
  colnames(eval_stats) <- c("auc", "kappa", "tss")
    
  pred_stats <- data.frame(matrix(ncol=6))
  colnames(pred_stats) <- c("n", "k", "ll", "AIC", "AICc", "BIC")
  
  predictions <- stack()
  
  varnum = length(names(env))
  variable_contributions <- data.frame(matrix(ncol = varnum))
  colnames(variable_contributions) <- c(names(env))
  
  count = 0
  
  for(i in 1:length(folders)){ #loop through model output folders
    
    directory <- paste(ModelOutputPath(),folders[i],sep="")
    setwd(directory)
    cat(sprintf('Folder: %s\n', folders[i]))

    for(j in 1:length(modelList[[i]]@models)-1){ #loop within folders among model outputs
      
      count = count + 1
      cat(sprintf('Replicate %s\n', j))
      
      presence <- read.csv(sprintf("species_%s_samplePredictions.csv", j))
      background <- read.csv(sprintf("species_%s_backgroundPredictions.csv", j))
      max_res<-read.csv("maxentResults.csv")
      lambdas <- sprintf("species_%s.lambdas", j)
      pp <- presence$Cloglog.prediction
      trainpp <- pp[presence$Test.or.train=="test"]
      bb <- background$Cloglog 

      #Variable Contribution/Permutation Importance
      # vnum = length(names(env))
      # conts <- vnum + 10 
      # imp <- conts + vnum
      # 
      # cat('Calculating variables')
      # contribution <- modelList[[i]]@results[,count][11:conts]
      # conts = conts + 1
      # importance <- modelList[[i]]@results[,count][conts:imp]
      # #how to solve this problem
      # #s0 <- modelList[[i]]@results[,i][index]
      # variables <- rbind(contribution, importance)
      # colnames(variables) <- names(env)
      # print(variables)
      
      #Model Prediction
      p <- predictMaxent(lambdas, env)
      predictions <- stack(predictions, p$prediction_cloglog)
      cat("\n")
      
      #Information Criteria
      infoc <- ic(p$prediction_raw, presences, lambdas)
      pred_stats[count,] <- infoc
      
      #AUC/ROC
      combined <- c(trainpp, bb)
      label <- c(rep(1,length(trainpp)), rep(0,length(bb)))
      pred <- prediction(combined,label)
      perf <- performance(pred, "tpr", "fpr")
      aucvalue <- performance(pred,"auc")@y.values[[1]]
      
      #AUC TEST / TRAIN / DIFF
      auctrain <- modelList[[i]]@results[,count][5]
      auctest <- modelList[[i]]@results[,count][8]
      aucsd <- modelList[[i]]@results[,count][9]
      aucdiff <- abs(auctest-auctrain)
      names(aucdiff) <- c("AUC.Difference")
      print(aucdiff)
      
      #Model Averaging Results
      # num_models <- ncol(modelList[[i]]@results)
      # average_contributions <- modelList[[i]]@results[,num_models][11:conts]
      # average_importance <- modelList[[i]]@results[,num_models][conts:imp]
      # average_variables <- rbind(average_contributions, average_importance)
      # print(average_variables)
      
      #KAPPA
      #Confusion matrix
      confusion <- function(thresh){
        return(cbind(c(length(trainpp[trainpp>=thresh]), length(trainpp[trainpp<thresh])),
                     c(length(bb[bb>=thresh]), length(bb[bb<thresh]))))
      }
      
      #Kappa statistic
      mykappa <- function(thresh){
        return(Kappa(confusion(thresh)))
      }
      
      #Calculate Kappa
      fpr = perf@x.values[[1]]
      tpr = perf@y.values[[1]]
      sum = tpr + (1-fpr)
      index = which.max(sum)
      cutoff = perf@alpha.values[[1]][[index]]
      kappavalue <- mykappa(cutoff)
      kappavalue <- mykappa(min(trainpp))
      kappa <- kappavalue$Unweighted[[1]]
      
      #TSS
      mytss <- function(thresh){
        con <- confusion(thresh)
        sensitivity<-con[1,1]/(con[1,1]+con[2,1])
        specificity<-con[2,2]/(con[2,2]+con[2,1])
        tss<-sensitivity+specificity-1
        return(tss)
      }
      
      #Calculate TSS
      tssvalue <- mytss(cutoff)
      
      eval_stats[count,] <- c(aucvalue, kappa, tssvalue)
      #cat(sprintf('AUC: %s, Kappa: %s, TSS: %s \n', aucvalue, kappa, tssvalue))

    }
    
  }
  return(list(eval_stats, pred_stats, predictions))
}

maxentVariableResults <- function(modelList){
  
  contribution <- maxent_assemble_results(modelList, include = c("auc", "contribution"))
  importance <- maxent_assemble_results(modelList, include = c("auc", "importance"))
  
  return(list(contribution, importance))
}

predictMaxent <- function(lambda, env){
  #Function for predicting maxent model outputs
  
  prediction <- rmaxent::project(lambda, env)
  
  return(prediction)
}

rasterAlgebra <- function(rasterStack){
  
  s <- stack(rasterStack)
  
  #raster difference
  #subtract
  
  #raster averaging
  rasterMean <- calc(s, mean)
  
  #raster median
  rasterMedian <- calc(s, median)
  
  #raster range
  rasterRange <- calc(s, range)
  
  #raster quantiles
  #Q25 <- calc(s, fun=stats::quantile, probs=.25, na.rm=TRUE)
  #Q75 <- calc(s, fun=stats::quantile, probs=.75, na.rm=TRUE)
  
  return(list(rasterMean, rasterMedian, rasterRange))
  
}

batchMaxent <- function(spp, env, iteration, cv){
  #Function to batch run maxent models for evaluation
  #Inputs: spp - species .csv file formatted as spp, lon, lat
  #######: env - climate .grd raster stack object
  #######: iteration - number of repeated model bootstraps
  #######: cv - number of kfolds for cross-validation
  
  pseudo <- sampleBackground(spp)
  PA <- extractClimateData(spp, env, pseudo)
  presence <- as.data.frame(PA[[1]])
  pseudoabsence <- as.data.frame(PA[[2]])
  
  CV <- crossValidation(presence, pseudoabsence, cv) 
  
  #initialize empty lists for model/prediction
  models <- list()
  count = 0

  for(i in 1:cv){ #number of cross validations
    
    train_p <- as.data.frame(CV[[i]][1:2])
    train_a <- as.data.frame(CV[[i]][3:4])
    test_p <- as.data.frame(CV[[i]][5:6])
    test_a <- as.data.frame(CV[[i]][7:8])
  
    for(j in 1:iteration){ #number of repeated model runs
      
      count = count + 1
      #model output directory
      opath = paste(ModelOutputPath(), sprintf("Model%s", count), sep="")
      cat(sprintf('Running Model Number %s \n', count))
      modelArgs <- maxentArgs()
      modelOut <- runMaxent(test_p, test_a, env, opath, modelArgs)
      models[[count]] <- modelOut
      
        if(count == cv*iteration){ 
          eval <- maxentStatistics(models, env, test_p) 
          cat('Finished evaluating models\n')
          print(cbind(eval[[1]], eval[[2]]))
          eval[[4]] <- rasterAlgebra(eval[[3]])
          #eval 1, 2 model statistics
          #eval 3 prediction stack
          #eval 4 algebra stack
          #model evaluation statistics
          #write.csv(eval, file = "../Model_Evaluation_Statistics.csv")
        }
    }
  }
return(eval)
}

#To Batch Run Maxent Models
MaxentBatch <- batchMaxent(subsamp.pa.all.coords.biomod, BioMod_reduced, iteration = 1, cv = 1)



###Code Testing



###Checklist

#remove cv from implementation? [x]
#implement GBM?
#run models for prism & climond data
#run model projections [x]
#calculate information criteria [x]
#calculate variable importance and permutation contribution [x]
#auc diff [x]
#raster output algebra [x]
#method to identify lowest scores/'best model'?
#maxent args function[x]
#add method to iterate through datasets
#add method to iterate through climate data


#R markdown
#complete r markdown script
#run 100 iterations of maxent models 
#add limiting function to 10 or so of these models
#calculate model range and mean statistics
#add github for batch maxent and gbm models













