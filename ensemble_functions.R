#regression errors 

#x estimate y true

rsq <- function(x,y) {
  1 - sum((y-x)^2)/sum((y-mean(y))^2)
}

#td test data , predicted test data
rmse <- function(td,ptd) {
  sqrt(mean((td - ptd)^2))
}



# combine the data and class for use with cluster functions, in order to keep the correct class (new)
combine_Data <- function(realdata,realclass){
  combined_data <- cbind(realdata,realclass)
  return(as.data.frame(combined_data))
  
} # combine_Data
#===

# shuffle the data and also export the shuffled class, in order to avoid overfitting (new)
shuffle_Data <- function(combine_Data, seed = 1234){
  set.seed(seed)
  shuffled_data <- sample_n(combine_Data,nrow(combine_Data))
  reslist <- list("shuffled_data" = shuffled_data, "shuffled_class" = shuffled_data[,ncol(shuffled_data)])
} # shuffle_data

# ==== cluster one time with dePDDP and export every step
onetime_clustering_export_instances <- function(normalized_red,adj = 0.5, kk = 40, tsamp_size = (nrow(normalized_red)*10)%/%100){
  
  #==================== cluster one time with idiv and export every step ===========================
  
  #=== clusterinf
  res<-dePDDP(normalized_red[,1:ncol(normalized_red)-1], k=kk, SP="DEN")
  
  #=== leafs occuring via final clustering
  leafs <- res[["leafs"]][["id"]]
  
  
  ## == get all node relations 
  allkids <- NULL
  kids <- NULL
  coln <- NULL
  for (i in 1 : length(res[["AllRes"]])) {
    if (res[["AllRes"]][[i]][["node"]][["Splitted"]] == 1) {
      kids[1] <- res[["AllRes"]][[i]][["node"]][["kid1.id"]]
      kids[2] <- res[["AllRes"]][[i]][["node"]][["kid2.id"]]
      
      
    }
    else{
      kids[1] <- 0
      kids[2] <- 0
    }
    allkids <- cbind(allkids,kids)
    coln <- c(coln,paste("node",i,sep = ""))
    colnames(allkids) <- coln
  }
  
  # ==
  
  
  # === create from relations the hypothetical instance for every number of clusters
  
  currentNodes <- NULL
  allNodes <- list()
  nodesInd <- NULL
  
  for (i in 1:length(leafs)) {
    if(i == 1) {
      
      currentNodes <- cbind(currentNodes,allkids[,1])
      allNodes[[1]] <- c(allNodes,currentNodes)
    }
    else{
      for (j in 1:length(allNodes[[i-1]])) {
        nodesInd <- c(nodesInd,allNodes[[i-1]][[j]])
      }  
      
      if(length(nodesInd[!(nodesInd %in% leafs)]) > 1) {
        namecol <- names(which.min(apply(allkids[,nodesInd[!(nodesInd %in% leafs)]],MARGIN=2,min)))
        replace_node <-  which( colnames(allkids)== namecol )
        nodesInd <- nodesInd[!(nodesInd == replace_node)]
        currentNodes <-  c(nodesInd, allkids[,replace_node])
      }
      else if (length(nodesInd[!(nodesInd %in% leafs)]) < 1) {
        print(paste("just ended in cluster number", i+1, sep = ""))
      }
      else if (length(nodesInd[!(nodesInd %in% leafs)]) == 1) {
        replace_node <- nodesInd[!(nodesInd %in% leafs)]
        nodesInd <- nodesInd[!(nodesInd == replace_node)]
        currentNodes <-  c(nodesInd, allkids[,replace_node])
      }
      
      nodesInd <- NULL
      allNodes [[i]] <- as.list(currentNodes)
      print(currentNodes)
    }
  }
  
  
  # === 
  
  datrow <- nrow(normalized_red)
  
  # ==== create the clustering ID matrix with the correct IDs (not repeated) for every nmber till the knum` that idiv returns
  clusterIDs_whole <- matrix(data = NA, nrow = datrow,ncol = length(allNodes))
  clusterIDs_whole <- as.data.frame(clusterIDs_whole)
  
  for (i in 1:length(allNodes)) {
    
    instance_clusters <- as.numeric(allNodes[[i]])
    print(instance_clusters)
    for (j in instance_clusters) {
      clusterIDs_whole[res[["AllRes"]][[j]][["node"]][["ids"]],i] <- j
    }
    
  }
  
  
  # ===== create the correct clid and test clid list
  cc <- list()
  for (i in 1:(length(allNodes)-1)) {
    ind <- i*2
    cc[[ind -1]] <- clusterIDs_whole[1:(datrow-tsamp_size),i]
    cc[[ind]] <- clusterIDs_whole[((datrow-tsamp_size)+1):datrow,i]
  }
  cc <- list("clustlist" = cc)
  # =====
  result <- list("cc" = cc, "leafs" = leafs)
  return(result)
}
#===



# Linear Regression with the schema, given the dataset and the clustering indexes  (new)
linRegEnsemble <- function(workingData_red, idExp, respVar, test_size = (nrow(workingData_red) * 10) %/% 100) {
  leafs <- idExp$leafs
  idExp <- idExp$cc
  nleafs <- length(leafs)
  
  datrow <- nrow(workingData_red)
  trainData <- workingData_red[1:(datrow - test_size),]
  # cc[[ind -1]] <- clusterIDs_whole[1:(datrow-tsamp_size),i]
  testdata <- workingData_red[((datrow - test_size) + 1):datrow,]
  testhealth <- testdata[, ncol(testdata)]
  
  #============= Linear regression rerun also
  
  estimates <- NULL
  rmses <- NULL
  rsqs <- NULL
  cnam <- NULL
  timeperclustnum <- NULL
  
  clusterlist <- idExp
  till_now <- NULL
  
  
  #check if response column has all the elements
  
  resCol <- respVar
  if (!resCol %in% colnames(workingData_red))
  {
    cat("Imput a correct column from the dataset!\n")
    
    return(0)
  }
  
  modelvars <- setdiff(colnames(workingData_red), resCol)
  fvar <-  paste(resCol, "~", sep = "")
  
  for (i in 2:nleafs) {
    #from already computed clusters
    ind <- i - 1
    show(i)
    #ksexwrise ta cluster samples gia na ftiakseis ta mondela
    if (ind == 1) {
      clusterID <- clusterlist[["clustlist"]][[ind]]
      clid <- table(clusterID)
      show(clid)
      testclusterID <- clusterlist[["clustlist"]][[ind + 1]]
    } else {
      clusterID <- clusterlist[["clustlist"]][[ind * 2 - 1]]
      clid <- table(clusterID)
      show(clid)
      testclusterID <- clusterlist[["clustlist"]][[ind * 2]]
    }
    
    cnam <- c(cnam, i)
    names <- NULL
    estimate <- NULL
    timeperrun <- NULL
    
    for (c in sort(unique(clusterID))) {
      if (!(c %in% till_now)) {
        TrainData <- trainData[which(clusterID[] == c),]
        
        rownames(TrainData) <- 1:nrow(TrainData)
        show(nrow(TrainData))
        nam <- paste("LM", c, sep = "")
        
        
        
        timeperclust <-
          system.time(assign(nam, lm(as.formula(
            paste(fvar, paste(modelvars, collapse = "+"))
          ), data = TrainData)))
        
        
        timeperrun <- cbind(timeperrun, timeperclust)
        till_now <- c(till_now, c)
        
      }
    }
    
    for (j in 1:nrow(testdata)) {
      nam2 <- paste("LM", testclusterID[j], sep = "")
      estimate[j] <-
        predict(eval(parse(text = nam2)), newdata = testdata[j,])
      
    }
    
    
    estimates <- cbind(estimates, estimate)
    colnames(estimates) <- cnam
    rmses <- cbind(rmses, rmse(estimate, testhealth))
    rsqs <-  cbind(rsqs, rsq(estimate, testhealth))
    colnames(rmses) <- cnam
    colnames(rsqs) <- cnam
    timeperclustnum <- cbind(timeperclustnum, sum(timeperrun[3,]))
    colnames(timeperclustnum) <- cnam
    
    gc()
  }
  
  ####plain estimation results (without ensemble yet)
  rmsesLRdppd_up <- rmses
  rsqsLRdppd_up <- rsqs
  
  timeperclustnumLRdppd_up <- timeperclustnum
  estimatesLRdppd_up <- estimates
  
  
  ####catholic LR
  TrainData <- workingData_red[1:(datrow - test_size),]
  nam <- paste("LR_all_3k")
  timeLRALL <-
    system.time(assign(nam, lm(as.formula(
      paste(fvar, paste(modelvars, collapse = "+"))
    ), data = TrainData)))
  
  
  allLRest <- NULL
  for (j in 1:nrow(testdata)) {
    allLRest[j] <-  predict(LR_all_3k, newdata = testdata[j,])
    
  }
  
  rmseCatholic <- rmse(allLRest, testhealth)
  rsqCatholic <- rsq(allLRest, testhealth)
  
  #### mean estimation LR (ensemble)
  name <- "estimatesLRdppd_up"
  
  #asc est,rmse
  stpoint <- 2
  cnam <- NULL
  testmeanest <- NULL
  testmeanrmseAcs <- NULL
  testmeanrsqAcs <- NULL
  for (finpoint in 3:length(leafs)) {
    cnam <- c(cnam, paste(stpoint, "to", finpoint, sep = ""))
    
    # est <- meanestimation_for_latestcheck(name, stpoint,finpoint)
    
    meanestimation <- NULL
    sp <- stpoint - 1
    fp <- finpoint - 1
    # show(fp)
    
    for (i in 1:test_size) {
      meanestimation[i] <- mean(estimatesLRdppd_up[i, sp:fp])
      
    }
    rmsemean <- rmse(meanestimation, testhealth)
    rsqmean <- rsq(meanestimation, testhealth)
    
    testmeanest <- cbind(testmeanest, meanestimation)
    testmeanrmseAcs <- cbind(testmeanrmseAcs, rmsemean)
    testmeanrsqAcs  <- cbind(testmeanrsqAcs, rsqmean)
    colnames(testmeanest) <- cnam
    colnames(testmeanrmseAcs) <- cnam
  }
  
  #######============lr dppd list updated===========
  lrlistdppd_up <- NULL
  tlr_up <- NULL
  tlr_up <- c(as.numeric(timeLRALL[3]), timeperclustnumLRdppd_up)
  tlr_up <- as.data.frame(tlr_up)
  tlr_up <- t(tlr_up)
  ch <- 2:length(leafs)
  colnames(tlr_up) <- c("timeLRALL", ch)
  
  
  lrlistdppd_up <-
    list(
      "times" = tlr_up,
      "estimates" = estimatesLRdppd_up,
      "rmses" = rmsesLRdppd_up,
      "rsqs" = rsqsLRdppd_up,
      "meanestimates_ascending" =  testmeanest,
      "rmses_mean_asc" = testmeanrmseAcs,
      "rsq_mean_asc" = testmeanrsqAcs
    )
  
  ## ++++++++++++++++++++++++++++++++++               timesums
  
  stpoint <- 2
  timeextra <- NULL
  cnamt <- NULL
  for (finpoint in 2:length(leafs)) {
    # show(finpoint)
    if (finpoint == 2) {
      ttime <- lrlistdppd_up[["times"]][finpoint]
    }
    else{
      ttime <-  lrlistdppd_up[["times"]][finpoint]
    }
    timeextra <- cbind(timeextra, ttime)
    
    cnamt <- c(cnamt, paste(stpoint, "to", finpoint, sep = ""))
    colnames(timeextra) <- cnamt
  }
  
  timesum <- NULL
  for (i in 2:length(leafs)) {
    # show(i)
    #cnamd <- c(cnamd,paste(stpoint,"to",finpoint,sep = ""))
    timesum <- cbind(timesum, sum(timeextra[1:i]))
  }
  
  
  timesum <- timesum[, 1:(length(leafs) - 2)]
  timesum <- t(as.data.frame(timesum))
  colnames(timesum) <- cnam
  timesumLRdppd_up_asc <- timesum
  
  
  estimatesLRdppd_up <- cbind(allLRest, estimatesLRdppd_up)
  rmsesLRdppd_up  <- cbind(rmseCatholic, rmsesLRdppd_up)
  rsqsLRdppd_up  <- cbind(rsqCatholic, rsqsLRdppd_up)
  
  ####update list
  lrlistdppd_up <-
    list(
      "times" = tlr_up,
      "estimates" = estimatesLRdppd_up,
      "rmses" = rmsesLRdppd_up,
      "rsqs" = rsqsLRdppd_up,
      "meanestimates_ascending" =  testmeanest,
      "rmses_mean_asc" = testmeanrmseAcs,
      "rsq_mean_asc" = testmeanrsqAcs,
      "time_mean_asc" = timesumLRdppd_up_asc
    )
  
  
  return(lrlistdppd_up)
  ####
  
} #linRegEnsemble
#===

# Random Forests regression with the schema, given the dataset and the clustering indexes  (new)
RFregEnsemble <- function(workingData_red, idExp, respVar, test_size = (nrow(workingData_red) * 10) %/% 100) {
  leafs <- idExp$leafs
  idExp <- idExp$cc
  nleafs <- length(leafs)
  
  datrow <- nrow(workingData_red)
  testdata <- workingData_red[((datrow - test_size) + 1):datrow, ]
  testhealth <- testdata[, ncol(testdata)]
  
  #### estimate their RF estimations
  
  estimates <- NULL
  rmses <- NULL
  rsqs <- NULL
  cnam <- NULL
  timeperclustnum <- NULL
  
  clusterlist <- idExp
  till_now <- NULL
  
  
  #check if response column has all the elements
  
  resCol <- respVar
  if (!resCol %in% colnames(workingData_red))
  {
    cat("Imput a correct column from the dataset!\n")
    
    return(0)
  }
  
  rFRNVars = length(colnames(workingData_red[1,])) - 1
  fvar <-  workingData_red[1:(datrow - test_size),resCol]
  
  #spliting the dataset into train and test according to the clusterIDs
  
  for (i in 2:nleafs) {
    #from already computed clusters
    ind <- i - 1
    show(i)
    #ksexwrise ta cluster samples gia na ftiakseis ta mondela
    if (ind == 1) {
      clusterID <- clusterlist[["clustlist"]][[ind]]
      clid <- table(clusterID)
      show(clid)
      testclusterID <- clusterlist[["clustlist"]][[ind + 1]]
    } else {
      clusterID <- clusterlist[["clustlist"]][[ind * 2 - 1]]
      clid <- table(clusterID)
      show(clid)
      testclusterID <- clusterlist[["clustlist"]][[ind * 2]]
    }
    
    #training of the models for each cluster
    
    cnam <- c(cnam, i)
    names <- NULL
    estimate <- NULL
    timeperrun <- NULL
    
    for (c in sort(unique(clusterID))) {
      if (!(c %in% till_now)) {
        rFRTrainData <- workingData_red[which(clusterID[] == c),]
        rownames(rFRTrainData) <- 1:nrow(rFRTrainData)
        show(nrow(rFRTrainData))
        # show(colnames(rFRTrainData))
        # show(anyNA(rFRTrainData))
        nam <- paste("RF_testsunolo", c, sep = "")
        #parrallel
        rFRNCores <- 6
        registerDoMC(rFRNCores) #number of cores on the machine
        tfvar <- rFRTrainData[,resCol]
        timeperclust <-
          system.time(rFRModel <-
                        foreach(y = seq(6), .combine = randomForest::combine) %dopar% {
                          rf <-
                            randomForest(tfvar ~ .,rFRTrainData,ntree = 5,mtry = ceiling(rFRNVars / 3),norm.votes = FALSE,importance = TRUE)
                        })
        timeperrun <- cbind(timeperrun, timeperclust)
        assign(nam, rFRModel)
        till_now <- c(till_now, c)
      }
    }
    
    #prediction of the dataset for a single clustering
    for (j in 1:nrow(testdata)) {
      nam2 <- paste("RF_testsunolo", testclusterID[j], sep = "")
      estimate[j] <-
        predict(eval(parse(text = nam2)), newdata = testdata[j,])
      
    }
    
    #collection of predictions for later tinkering
    estimates <- cbind(estimates, estimate)
    colnames(estimates) <- cnam
    rmses <- cbind(rmses, rmse(estimate, testhealth))
    rsqs <-  cbind(rsqs, rsq(estimate, testhealth))
    colnames(rmses) <- cnam
    colnames(rsqs) <- cnam
    timeperclustnum <- cbind(timeperclustnum, sum(timeperrun[3,]))
    colnames(timeperclustnum) <- cnam
    
    gc()
  }
  
  ####plain estimation results (without ensemble yet)
  rmsesRFdpddp_up <- rmses
  rsqsRFdpddp_up <- rsqs
  timeperclustnumRFdpddp_up <- timeperclustnum
  estimatesRFdpddp_up <- estimates
  
  
  ####catholic RF training
  rFRNCores <- 6
  registerDoMC(rFRNCores) #number of cores on the machine
  timeRFALL <-
    system.time(rFRModel <-
                  foreach(y = seq(6), .combine = randomForest::combine) %dopar% {
                    rf <-
                      randomForest(
                        fvar ~ .,
                        workingData_red[1:(datrow - test_size),],
                        ntree = 5,
                        mtry = ceiling(rFRNVars / 3),
                        norm.votes = FALSE,
                        importance = TRUE
                      )
                  })
  
  
  assign("RF_all_3k", rFRModel)
  
  ####prediction of the dataset for the catholic model (maybe slow way, though not of significance for this study)
  allRFest <- NULL
  for (j in 1:nrow(testdata)) {
    allRFest[j] <-  predict(rFRModel, newdata = testdata[j,])
    
  }
  
  rmseCatholic <- rmse(allRFest, testhealth)
  rsqCatholic <- rsq(allRFest, testhealth)
  
  #### mean estimation RF for regression (ensemble) 
  
  #asc est,rmse
  stpoint <- 2
  cnam <- NULL
  testmeanest <- NULL
  testmeanrmseAcs <- NULL
  testmeanrsqAcs <- NULL
  for (finpoint in 3:length(leafs)) {
    cnam <- c(cnam, paste(stpoint, "to", finpoint, sep = ""))
    
    meanestimation <- NULL
    sp <- stpoint - 1
    fp <- finpoint - 1
    # show(fp)
    
    for (i in 1:test_size) {
      meanestimation[i] <- mean(estimatesRFdpddp_up[i, sp:fp])
      
    }
    rmsemean <- rmse(meanestimation, testhealth)
    rsqmean <- rsq(meanestimation, testhealth)
    
    testmeanest <- cbind(testmeanest, meanestimation)
    testmeanrmseAcs <- cbind(testmeanrmseAcs, rmsemean)
    testmeanrsqAcs  <- cbind(testmeanrsqAcs, rsqmean)
    colnames(testmeanest) <- cnam
    colnames(testmeanrmseAcs) <- cnam
  }
  
  
  
  ####============rf dppd list updated===========
  rflistdppd_up <- NULL
  
  #time colection
  trf_up <- NULL
  trf_up <- c(as.numeric(timeRFALL[3]), timeperclustnumRFdpddp_up)
  trf_up <- as.data.frame(trf_up)
  trf_up <- t(trf_up)
  ch <- 2:length(leafs)
  colnames(trf_up) <- c("timeRFALL", ch)
  
  rflistdppd_up <-
    list(
      "times" = trf_up,
      "estimates" = estimatesRFdpddp_up,
      "rmses" = rmsesRFdpddp_up,
      "rsqs" = rsqsRFdpddp_up,
      "meanestimates_ascending" =  testmeanest,
      "rmses_mean_asc" = testmeanrmseAcs,
      "rsq_mean_asc" = testmeanrsqAcs
    )
  
  
  ## ++++++++++++++++++++++++++++++++++               timesums
  
  
  stpoint <- 2
  timeextra <- NULL
  cnamt <- NULL
  for (finpoint in 2:length(leafs)) {
    # show(finpoint)
    if (finpoint == 2) {
      ttime <- rflistdppd_up[["times"]][finpoint]
    }
    else{
      ttime <-  rflistdppd_up[["times"]][finpoint]
    }
    timeextra <- cbind(timeextra, ttime)
    
    cnamt <- c(cnamt, paste(stpoint, "to", finpoint, sep = ""))
    colnames(timeextra) <- cnamt
  }
  
  timesum <- NULL
  for (i in 2:length(leafs)) {
    # show(i)
    timesum <- cbind(timesum, sum(timeextra[1:i]))
  }
  
  timesum <- timesum[, 1:(length(leafs) - 2)]
  timesum <- t(as.data.frame(timesum))
  colnames(timesum) <- cnam
  
  timesumRFdppd_up_asc <- timesum
  estimatesRFdpddp_up <- cbind(allRFest, estimatesRFdpddp_up)
  rmsesRFdpddp_up  <- cbind(rmseCatholic, rmsesRFdpddp_up)
  rsqsRFdpddp_up  <- cbind(rsqCatholic, rsqsRFdpddp_up)
  
  
  ####update list
  rflistdppd_up <-
    list(
      "times" = trf_up,
      "estimates" = estimatesRFdpddp_up,
      "rmses" = rmsesRFdpddp_up,
      "rsqs" = rsqsRFdpddp_up,
      "meanestimates_ascending" =  testmeanest,
      "rmses_mean_asc" = testmeanrmseAcs,
      "rsq_mean_asc" = testmeanrsqAcs,
      "time_mean_asc" = timesumRFdppd_up_asc
    )
  
  return(rflistdppd_up)
  ####
  
} #RFregEnsemble
#===