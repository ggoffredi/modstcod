#Multivariate Spatio Temporal Weighted Anomaly Detector, with Linear Fisher Subset Scan


######## Main Algorithm ##########################################

#cellData: list of n kpi, each containing a matrix for all element
#neigh_number: includes self in count of neighborhood
MuST <- function(cellData, trainingData = NULL, coordinates, festive, period = 24, neigh_number = 20, visual = FALSE, 
                 methodList = NULL, anomaly_window = 4, NA_strat = "local_mean",
                 method_weighting = TRUE, best = TRUE, feature_set = FALSE) {
  
  state_mem <- 5
  
  memory <- 24
  
  confidence_threshold <- 0.997
    
  min_difference <- 0.0001
  
  kpi_number <- length(cellData)
  node_number <- nrow(cellData[[1]])
  final_timeslot <- ncol(cellData[[1]])
  
  
  
  diff_lag <- 3
  diff_differences <- 1
  diff_threshold <- 3
  diff_thr_conf <- 0.68
  
  min_diff_window <- (diff_lag * diff_differences) + 1
  
  
  #method info
  
  if(is.null(methodList)) {
    method_number <- 6
  }
  else {
    method_number  <- length(methodList$active)
  }
  selectMethod <- vector("logical", length = method_number)
  startMethod <- vector("numeric", length = method_number)
  nameMethod <- vector("character", length = method_number)
  
  
  if(is.null(methodList)) {
    selectMethod <- rep(TRUE, method_number)
    nameMethod <- c("season", "dist", "trend","last", "history", "cdf", "diff", "diff_season")
    startMethod <- c(1,1,2,2,1,1, min_diff_window, memory + 1)
  }
  else {
    selectMethod <- methodList$active
    nameMethod <- methodList$names
    startMethod <- methodList$starts
  }
  
  
  threshold_confidence <- 0.68
  
  #indicates quality of method
  method_weight <- rep(list(matrix(1, nrow = method_number, ncol = kpi_number)), node_number)
  
  
  
  anomaly_significance <- 0.997
  max_priority <- 0.999999
  
  
  
  #following data is stored in list instead of matrix to simulate cell distributed behavior
  #element in a list[[1]] represents knowledge of node 1
  #this is redundant in this algorithm, as cluster_season_sample should be the same
  
  #seasonality data
  cell_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                      nrow = kpi_number, ncol = period)), node_number)
  
  cell_F_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                        nrow = kpi_number, ncol = period)), node_number)
  
  cell_seasonality_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                         nrow = kpi_number, ncol = period)), node_number)
  cell_F_seasonality_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                           nrow = kpi_number, ncol = period)), node_number)
  
  
  #previous values
  cell_memories <- rep(list(matrix(rep(0, kpi_number * memory), 
                                   nrow = kpi_number, 
                                   ncol = memory)), node_number)
  
  #historical values
  cell_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  cell_diff_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_diff_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  cell_diff_periodical_mean <- rep(list(matrix(rep(0, kpi_number*period), 
                                               nrow = kpi_number, ncol = period)), node_number)
  cell_diff_periodical_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                             nrow = kpi_number, ncol = period)), node_number)
  
  
  cell_cdf <- rep(list(rep(list(NULL), kpi_number)), node_number)
  
  
  cell_descent_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_descent_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_distance_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_distance_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  neigh_index <- rep(list(vector("numeric", length = neigh_number)), node_number)
  
  prev_neigh_value <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
  prev_priority <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
  neigh_weight <- rep(list(vector("numeric", length = neigh_number)), node_number)
  
  node_votes <- vector("numeric", length = node_number)
  final_result <- NULL
  
  #priority in neighborhood
  #node_priority[[node]][1] refers to same node
  #node_priority[[node]][2:neigh_number] refers to neighbours of node
  #righ now no order between neighbours is to be assumed
  node_priority <- rep(list(vector("numeric", length = neigh_number)), node_number)
  
  anomaly_collection <- rep(list(rep(list(NULL), node_number)), kpi_number)
  
  anomaly_collection2 <- rep(list(NULL),  node_number)
  
  anomaly_times <- rep(list(NULL), node_number)
  
  #because this implementation is not distributed, to avoid computing n*k times a bestK clustering
  #only the first node does, and then it's assumed valid for all
  
  bestK <- NULL
  bestLevels <- rep(list(NULL), node_number)
  
  #################### ALGORITHM
  
  skipTraining <- FALSE
  
  if(!is.null(trainingData)) {
    skipTraining <- TRUE
    
    cell_seasonality <- trainingData$seasonality
    cell_F_seasonality <- trainingData$F_seasonality
    cell_seasonality_sd <- trainingData$season_sd
    cell_F_seasonality_sd <- trainingData$F_season_sd
    cell_historical_mean <- trainingData$hist_mean
    cell_historical_sd <- trainingData$hist_sd
    cell_F_historical_mean <- trainingData$F_hist_mean
    cell_F_historical_sd <- trainingData$F_hist_sd
    cell_cdf <- trainingData$kernel_cdf
    
    if(method_weighting == TRUE) {
      method_weight <- trainingData$method_weights
    }
    cell_descent_historical_mean <- trainingData$descent_mean
    cell_descent_historical_sd <- trainingData$descent_sd
    cell_distance_historical_mean <- trainingData$distance_mean
    cell_distance_historical_sd <- trainingData$distance_sd
    cell_diff_mean <- trainingData$diff_mean
    cell_diff_sd <- trainingData$diff_sd
    cell_diff_periodical_mean <- trainingData$diff_season_mean
    cell_diff_periodical_sd <- trainingData$diff_season_sd
    
    
  }
  
  
  
  #NA_STRATEGY
  if(NA_strat == "zero") {
    
    #for simplicity test without na
    for(feature in 1:length(cellData)) {
      cellData[[feature]][is.na(cellData[[feature]])] <- 0
      
    }
  }
  else if (NA_strat == "local_mean") {
    for(feature in 1:length(cellData)) {
      
      cellData[[feature]] <- halfNA(cellData[[feature]])
      
      #what can't be imputed put to 0
      cellData[[feature]][is.na(cellData[[feature]])] <- 0
    }    
    
  }
  
  #stop using data.frame
  for(feature in 1:length(cellData)) {
    cellData[[feature]] <- as.matrix(cellData[[feature]])
  }
  
  #initialize neighbors
  
  
  for(node in (1:node_number)) {
    
    library(FNN)
    #should use matrix of spatial coord, ordered by node number
    #as a test, just use the index of the node    
    if(is.null(coordinates)) {
      coordinates <- cbind(c(1:node_number), c(1:node_number))
    }
    spnn <- get.knn(coordinates, neigh_number - 1, "kd_tree")
    nn.index <- spnn[[1]]
    
    neigh_index[[node]] <- c(node, nn.index[node,])
    
    #For starter, initialize neigh_weight to 1
    neigh_weight[[node]] <- rep(1, length(neigh_weight[[node]]))
    
    
  }
  
  
  anomaly_clusters <- rep(list(NULL), final_timeslot)
  
  
  
  for(time in (1:final_timeslot)) {
    
    print(paste("Time: ", time))
    
    node_votes <- rep(0, node_number)
    
    #phase 1 - read and self detection
    
    ## START NODE OPERATION ##
    
    for(node in (1:node_number)) {
      
      
      #check time slot
      curr_period <- ((time-1) %% period) + 1
      
      #check if festive
      curr_day <- floor((time-1)/period) + 1
      festivity <- festive[curr_day]
      
      
      curr_values <- vector("numeric", length = kpi_number)
      
      method_fit <- rep(list(vector("numeric", length = kpi_number)) ,method_number)
      
      
      
      ########################Parameter preparation
      
      #for each feature in the node, calculate various measures of unlikelihood
      for(kpi in (1:kpi_number)) {
        
        curr_values[kpi] <- cellData[[kpi]][node, time]
        
        if(festivity) {
          curr_seasonality <- cell_F_seasonality[[node]][kpi,]
          
        }
        else {
          curr_seasonality <- cell_seasonality[[node]][kpi,]
          
        }
        
        
        #find standardized values for distance and trend
        season_avg <- mean(curr_seasonality)
        season_rng <- (range(curr_seasonality)[2] 
                       - range(curr_seasonality)[1]) / 2
        
        if(season_rng < min_difference ) {
          season_rng <- min_difference
        }
        
        curr_val_std <- (curr_values[kpi] - season_avg) / season_rng 
        curr_season_std <- (curr_seasonality[curr_period] - season_avg) / season_rng 
        
        dist <- curr_val_std - curr_season_std
        
        if(curr_period > 1) {
          last_period <- curr_period - 1 
        }
        else {
          last_period <- period
        }          
        
        
        if(time >= 2) {
          
          
          prev_val_std <- (cell_memories[[node]][kpi,ncol(cell_memories[[node]])] - season_avg) / season_rng 
          prev_season_std <- (curr_seasonality[last_period] - season_avg) / season_rng 
          
          residual <- (curr_val_std - prev_val_std) - (curr_season_std - prev_season_std)
          
        }
        
        if(time >= min_diff_window)  {
          
          last <- ncol(cell_memories[[node]])
          first <- last - min_diff_window + 2
          diff_memories <- c(cell_memories[[node]][kpi, first:last], curr_values[kpi])
          
          last <- curr_period
          first <- last - min_diff_window + 1
          if(first < 1) {
            
            window2 <- curr_seasonality[1:last]
            #check if previous window is festive or not
            first <- first + 24
            if(festive[curr_day -1]) {
              window1 <- cell_F_seasonality[[node]][kpi, first:period]
            }
            else {
              window1 <- cell_seasonality[[node]][kpi, first:period]
            }
            
            diff_season <- c(window1, window2)
          }
          else {
            diff_season <- curr_seasonality[first:last]
          }
          
          diff_memories <- (diff_memories - season_avg) / season_rng
          diff_season <- (diff_season - season_avg) / season_rng
          
        }
        
        if(cell_historical_sd[[node]][kpi] < min_difference) {
          cell_historical_sd[[node]][kpi] <- min_difference
        }
        
        
        ######################### Method application
        
        for (meth in 1:method_number) {
          
          if(time >= startMethod[meth] && selectMethod[meth]) {
            
            if(nameMethod[meth] == "season") {
              
              method_fit[[meth]][kpi] <- featureSeasonFitness(curr_seasonality,
                                                              curr_values[kpi],
                                                              curr_period, replicas = 3, visual = visual)
              
            }            
            else if(nameMethod[meth] == "dist" ) {
              
              method_fit[[meth]][kpi] <- featureDistanceThresholdPVal(curr_season_std,
                                                                      curr_val_std,
                                                                      cell_distance_historical_sd[[node]][kpi],
                                                                      threshold_confidence, visual)
              
            }
            else if(nameMethod[meth] == "trend") {
              method_fit[[meth]][kpi] <- featureTrendThresholdPVal(curr_season_std,
                                                                   prev_season_std,
                                                                   curr_val_std,
                                                                   prev_val_std,
                                                                   cell_descent_historical_sd[[node]][kpi],
                                                                   threshold_confidence, visual) 
            }
            else if(nameMethod[meth] == "last") {
              
              method_fit[[meth]][kpi] <- featureLastValuesFitness(cell_memories[[node]][kpi,],
                                                                  curr_values[kpi], visual)
            }
            else if(nameMethod[meth] == "history") {
              
              if(festivity) {
                method_fit[[meth]][kpi] <- featureHistoryFitness(cell_F_historical_mean[[node]][kpi],
                                                                 cell_F_historical_sd[[node]][kpi],
                                                                 curr_values[kpi], visual)
                
              }
              else {
                method_fit[[meth]][kpi] <- featureHistoryFitness(cell_historical_mean[[node]][kpi],
                                                                 cell_historical_sd[[node]][kpi],
                                                                 curr_values[kpi], visual)
                
              }
              
            }
            else if(nameMethod[meth] == "cdf") {
              
              method_fit[[meth]][kpi] <- featureResidualCDFPVal(curr_seasonality[curr_period] - curr_values[kpi],
                                                                cell_cdf[[node]][[kpi]])
              
              
            }
            
            else if(nameMethod[meth] == "diff") {
              #print(nameMethod[meth])
              method_fit[[meth]][kpi] <- featureDiffPVal(diff_season,
                                                         diff_memories,
                                                         cell_diff_sd[[node]][kpi],
                                                         diff_thr_conf,
                                                         diff_lag,
                                                         diff_differences)
              
            }            
            else if(nameMethod[meth] == "diff_season") {
              #print(nameMethod[meth])
              method_fit[[meth]][kpi] <- featureDiffPVal(diff_season,
                                                         diff_memories,
                                                         cell_diff_periodical_sd[[node]][kpi,curr_period],
                                                         diff_thr_conf,
                                                         diff_lag,
                                                         diff_differences)
              
            } 
            
            
          }
          
          
        }
        
        
      }
      
      ############## Combining methods
      
      #CHECK THIS
      base <- 0.6
      
      #aggregate by kpi
      univariateEnsemblePriority <- vector("numeric", length = kpi_number)
      
      for(kpi in (1:kpi_number)) {
        
        kpiPVals <- NULL
        for(meth in 1:method_number) {
          
          if(time >= startMethod[meth] && selectMethod[meth]) {
            
            wei_po <- 1 - priorityWeight(1 - method_fit[[meth]][kpi],
                                         method_weight[[node]][meth,kpi],
                                         base)
            #weigh methods by quality
            kpiPVals <- c(kpiPVals,
                          wei_po)
            
          }
          
          
          
        }
        
        univariateEnsemblePriority[kpi] <- fisherCombinedPriority(kpiPVals)
        
      }
      
      #aggregate by method
      multivariateMethodPriority <- NULL 
      
      for(meth in 1:method_number) {
        
        if(time >= startMethod[meth] && selectMethod[meth]) {
          multivariateMethodPriority <- c(fisherCombinedPriority(method_fit[[meth]]))
          
        }
        
        
        
      }
      
      if(any(multivariateMethodPriority > confidence_threshold)) {
        
        for( method in (1: length(multivariateMethodPriority))) {
          
          if(multivariateMethodPriority[method] > confidence_threshold) {
            
            #print(paste("Anomaly node ",node," time ", time," signaled by method ", method ," over all kpi"))
            
          }
          
        }
        
        
      }
      else if(any(univariateEnsemblePriority > confidence_threshold)) {
        #print(paste("Anomaly node ",node," time ", time," signaled over a kpi by ensemble of methods"))
        
      }
      alfa <- confidence_threshold
      
      node_priority[[node]][1] <- priority_function(1 - univariateEnsemblePriority, alfa)
      
      prev_priority[[node]][1,] <-  shiftRight(prev_priority[[node]][1,] , 1)
      prev_priority[[node]][1,1] <- node_priority[[node]][1]
      
      if(node_priority[[node]][1] > anomaly_significance) {
        
        if(feature_set == TRUE) {
          
          kpis <- paste(which(univariateEnsemblePriority > confidence_threshold), collapse = ", ")
          print(paste("Anomalous node ",node,", kpis ",kpis))
        }
        else {
          print(paste("Anomalous node ",node,", overall"))
          
        }
      }
      
      
      #remember values
      cell_memories[[node]] <- cbind(cell_memories[[node]][,2:ncol(cell_memories[[node]])], curr_values)
      
      
    }
    
    ## END NODE OPERATION ##
    
    
    
    
    #Before cluster operation, synch, all nodes finish phase 1
    
    curr_period <- ((time-1) %% period) + 1
    
    
    
    #Propagation phase: each node sends its neighbors within a given range its current status
    #and value    
    
    if(neigh_number > 1) {
      
      for(node in (1:node_number)) {
        for(neigh_ind in (2:neigh_number)) {
          
          #index to read from main table
          true_neighbour_index <- neigh_index[[node]][neigh_ind]
          
          prev_priority[[node]][neigh_ind,] <- shiftRight(prev_priority[[node]][neigh_ind,], 1)
          prev_priority[[node]][neigh_ind,1] <- prev_priority[[true_neighbour_index]][1,1]
          
          node_priority[[node]][neigh_ind] <- node_priority[[true_neighbour_index]][1]
          
        }
      }
    }
  
    #START CLUSTER OPERATION
    
    
    for(node in (1:node_number)) {
      
      curr_period <- ((time-1) %% period) + 1
      
      
      #VOTING PROCEDURE
      
      if(time >= max(1 + 2 + ncol(prev_priority[[node]]))) {
        
        
        cluster_op <- fastGeneralWeightedFisherScan(node_priority[[node]], neigh_weight[[node]], 
                                                    anomaly_significance, best = best)
        
        if (cluster_op$vote > anomaly_significance) {
          
          #print(paste("Node: ",node,", anomaly cluster at time ",time))
          
          
          inv_neighs <- cluster_op$nodes
          inv_nodes <- NULL
          
          if(length(inv_neighs) >= 2) {
            
            
            for(anom_neigh  in inv_neighs) {
              inv_nodes <- c(inv_nodes, neigh_index[[node]][anom_neigh])
              
            }
            
            #Qui c'era il passaggio di node_votes
            
            node_votes[inv_nodes] <- node_votes[inv_nodes] + 1
            
            string <- paste(inv_nodes, collapse = " ")
            print(paste("(",node,",",time,")Involved nodes: ", string,", vote ",cluster_op$vote))
            
            new_ac <- anomaly_cluster.merge(anomaly_clusters[[time]], inv_nodes, cluster_op$vote)
            
            if(!is.null(new_ac)) {
              anomaly_clusters[[time]] <- new_ac
            }
            
            #only record anomaly if includes this node
            if(node %in% inv_nodes) {
              
              
              anomaly_times[[node]] <- c(anomaly_times[[node]], time)
              
            }
            
          }
          
        }
        
      }
      
      
      #SIMILARITY CHECK - WEIGHT UPDATING
      
      if(time >= 2 + ncol(prev_priority[[node]])) {
        similarity_update_rate <- 0.2
        
        for(neigh_ind in (2:neigh_number)) {
          
          
          similarity <- similarity_check(prev_priority[[node]][1,], prev_priority[[node]][neigh_ind,])
          
          neigh_weight[[node]][neigh_ind] <- ((1 - similarity_update_rate)*neigh_weight[[node]][neigh_ind] 
                                              + similarity_update_rate * similarity)
          
        }
        
        
      }
      
      
      
    }
    #END CLUSTER OPERATION
    
    #check votes
    if(!is.null(anomaly_clusters[[time]])) {
      removeX <- NULL
      for(clustIndex in 1:length(anomaly_clusters[[time]])) {
        clust <- anomaly_clusters[[time]][[clustIndex]]
        inv_nodes <- clust$nodes
        meaningful_nodes <- inv_nodes[which(node_votes[inv_nodes] > 2)]
        
        if(length(meaningful_nodes) <= 0) {
          
          removeX <- c(removeX, clustIndex)
          
        }
        anomaly_clusters[[time]][[clustIndex]]$nodes <- meaningful_nodes
      }
      
      for(remInd in rev(removeX)) {
        
        anomaly_clusters[[time]][[remInd]] <- NULL
        
      }
      if(length(anomaly_clusters[[time]]) <= 0) {
        anomaly_clusters[time] <- list(NULL)
      }
      
      for(clust in anomaly_clusters[[time]]) {
        string <- paste(clust$nodes, collapse = ", ")
        print(paste("Cluster nodes: ",string,"; priority ",clust$priority))
      }
    }
    
    
    final_result <- cbind(final_result, node_votes)
    
  }
  
  print("end")
  
  
  
  
  ################### Anomaly collection
  
  for(node in 1:node_number) {
    
    anomaly_times <- which(final_result[node,] > 0)
    
    motif_set <- NULL
    current_motif <- NULL
    
    for(cand in anomaly_times) {
      #create vector with data in anomaly time
      for(feat in cellData) {
        current_motif <- cbind(current_motif, feat[node,cand])
      }
      
      #add lines before and after anomaly time
      before <- 0
      after <- 0
      noMore <- F
      init_length <- nrow(current_motif)
      
      min_window <- anomaly_window
      
      while(nrow(current_motif) < min_window && noMore == FALSE) {
        #pad data
        
        if(before <= after && (cand - init_length - before) >= 1 ) {
          
          #add line before current
          
          line <- NULL
          for(feat in cellData) {
            line <- cbind(line,feat[node, cand - init_length - before])
            
          }
          
          current_motif <- rbind(line, current_motif)
          
          before <- before + 1
        }
        else if(cand + after + 1 <= ncol(cellData[[1]]))   {
          
          #add line after current
          
          line <- NULL
          for(feat in cellData) {
            line <- cbind(line,feat[node, cand + after + 1])
            
          }
          
          current_motif <- rbind(current_motif, line)
          
          after <- after + 1
        }
        else if((cand - init_length - before) >= 1 ) {   
          
          #add line before current
          
          line <- NULL
          for(feat in cellData) {
            line <- cbind(line,feat[node, cand - init_length - before])
            
          }
          
          current_motif <- rbind(line, current_motif)
          
          before <- before + 1
        }
        else {
          noMore <- TRUE
        }
      }
      
      motif_set <- c(motif_set, list(current_motif))
      current_motif <- NULL
      
    }
    
    
    if(!is.null(motif_set)) {
    anomaly_collection2[[node]] <- motif_set
    }
    
#     for(feature in 1:length(cellData)) {
#       found_anomalies <- anomaly_window_search(as.numeric(cellData[[feature]][node,]), anomaly_times, anomaly_window)
#       if(!is.null(found_anomalies)) {
#         anomaly_collection[[feature]][[node]] <- found_anomalies
#       }
#     }
    
  }
  
  return(list(result = final_result,
              anomalies = anomaly_collection2,
              clusters = anomaly_clusters))
  
}




######## Training Algorithm ##############



MuST.training <- function(cellData, festive, coordinates, period = 24, neigh_number = 10,
                          NA_strat = "local_mean",
                          season = FALSE, 
                          calc_cdf = FALSE,
                          fits = FALSE,
                          weighting = FALSE,
                          preCalcs = NULL) {
  #for kernel cdf
  library(sROC)
  
  state_mem <- 5
  
  confidence_threshold <- 0.997
  
  memory <- 7
  
  min_difference <- 0.0001
  
  kpi_number <- length(cellData)
  node_number <- nrow(cellData[[1]])
  final_timeslot <- ncol(cellData[[1]])
  
  diff_lag <- 3
  diff_differences <- 1
  diff_thr_conf <- 0.68
  
  min_diff_window <- (diff_lag * diff_differences) + 1
  
  
  
  nameMethod <- c("trend","cdf","diff","diff_season", "last")
  method_number = length(nameMethod)
  selectMethod <- rep(TRUE, method_number)
  startMethod <- c(2,1,min_diff_window,min_diff_window, memory + 1)  
  
  threshold_confidence <- 0.842
  
  #indicates quality of method
  method_weight <- rep(list(matrix(1, nrow = method_number, ncol = kpi_number)), node_number)
  
  methWeightPLis <-  rep(list(rep(list(rep(list(vector("numeric", length = 0)), method_number)), kpi_number)), node_number)
  methFits <-  rep(list(rep(list(rep(list(vector("numeric", length = 0)), method_number)), kpi_number)), node_number)
  
  anomaly_significance <- 0.9545
  max_priority <- 0.999999
  
  
  
  #following data is stored in list instead of matrix to simulate cell distributed behavior
  #element in a list[[1]] represents knowledge of node 1
  #this is redundant in this algorithm, as cluster_season_sample should be the same
  
  #seasonality data
  cell_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                      nrow = kpi_number, ncol = period)), node_number)
  cell_F_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                        nrow = kpi_number, ncol = period)), node_number)
  
  cell_seasonality_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                         nrow = kpi_number, ncol = period)), node_number)
  cell_F_seasonality_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                           nrow = kpi_number, ncol = period)), node_number)
  
  cell_STL_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                          nrow = kpi_number, ncol = period)), node_number)
  cell_F_STL_seasonality <- rep(list(matrix(rep(0, kpi_number*period), 
                                            nrow = kpi_number, ncol = period)), node_number)
  
  
  
  #previous values
  cell_memories <- rep(list(matrix(rep(0, kpi_number *memory), 
                                   nrow = kpi_number, 
                                   ncol = memory)), node_number)
  
  
  
  cell_cdf <- rep(list(rep(list(NULL), kpi_number)), node_number)
  
  
  
  #historical values
  cell_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_F_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  cell_descent_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_descent_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_diff_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_diff_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_diff_periodical_mean <- rep(list(matrix(rep(0, kpi_number*period), 
                                               nrow = kpi_number, ncol = period)), node_number)
  cell_diff_periodical_sd <- rep(list(matrix(rep(0, kpi_number*period), 
                                             nrow = kpi_number, ncol = period)), node_number)
  
  cell_distance_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_distance_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  neigh_index <- rep(list(vector("numeric", length = neigh_number)), node_number)
  
  prev_neigh_value <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
  prev_priority <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
  neigh_weight <- rep(list(vector("numeric", length = neigh_number)), node_number)
  
  node_votes <- vector("numeric", length = node_number)
  final_result <- NULL
  
  #priority in neighborhood
  #node_priority[[node]][1] refers to same node
  #node_priority[[node]][2:neigh_number] refers to neighbours of node
  #righ now no order between neighbours is to be assumed
  node_priority <- rep(list(vector("numeric", length = neigh_number)), node_number)
  anomaly_collection <- rep(list(NULL), node_number)
  anomaly_times <- rep(list(NULL), node_number)
  
  #because this implementation is not distributed, to avoid computing n*k times a bestK clustering
  #only the first node does, and then it's assumed valid for all
  
  bestK <- NULL
  bestLevels <- rep(list(NULL), node_number)
  
  
  ###################Algorithm start
  
  
  
  print("Imputing NA...")
  #NA_STRATEGY
  if(NA_strat == "zero") {
    
    #for simplicity test without na
    for(feature in 1:length(cellData)) {
      cellData[[feature]][is.na(cellData[[feature]])] <- 0
    }
  }
  else if (NA_strat == "local_mean") {
    for(feature in 1:length(cellData)) {
      
      cellData[[feature]] <- halfNA(cellData[[feature]])
      cellData[[feature]] <- halfNaN(cellData[[feature]])        
      
    }    
    
  }
  
  #stop using data.frame
  for(feature in 1:length(cellData)) {
    cellData[[feature]] <- as.matrix(cellData[[feature]])
  }
  
#   #initialize neighbors
#   print("Neighbor initialization...")
#   for(node in (1:node_number)) {
#     
#     library(FNN)
#     #should use matrix of spatial coord, ordered by node number
#     #as a test, just use the index of the node    
#     if(is.null(coordinates)) {
#       coordinates <- cbind(c(1:node_number), c(1:node_number))
#     }
#     spnn <- get.knn(coordinates, neigh_number - 1, "kd_tree")
#     nn.index <- spnn[[1]]
#     
#     neigh_index[[node]] <- c(node, nn.index[node,])
#     
#     #For starter, initialize neigh_weight to 1
#     neigh_weight[[node]] <- rep(1, length(neigh_weight[[node]]))
#     
#     
#   }
#   
  
  
  #calculate mean and sd for every node
  
  print("Computing seasonalities...")
  
  for(node in (1:node_number)) {
    for(kpi in (1:kpi_number)) {
      
      data <- cellData[[kpi]][node, ]
      des_data <- cellData[[kpi]][node,]
      
      
      diffy <- diff(data, diff_lag, diff_differences)
      #CHECK
      diffy <- c(rep(mean(diffy), diff_lag * diff_differences), diffy)   
      
      if(season == TRUE) {
        
        #seasonality
        for(hour in (0:(period-1))) {
          svec <- seq(1 + hour, length(data), period)
          festive_ind <- which(festive == T)
          
          cell_F_seasonality[[node]][kpi,hour + 1] <- mean(data[svec[festive_ind]], na.rm = T)
          cell_F_seasonality_sd[[node]][kpi,hour + 1] <- sd(data[svec[festive_ind]], na.rm = T)
          
          cell_seasonality[[node]][kpi,hour + 1] <- mean(data[svec[-festive_ind]], na.rm = T)
          cell_seasonality_sd[[node]][kpi,hour + 1] <- sd(data[svec[-festive_ind]], na.rm = T)
          
          cell_diff_periodical_mean[[node]][kpi,hour + 1] <- mean(diffy[svec], na.rm = T)
          cell_diff_periodical_sd[[node]][kpi,hour + 1] <- sd(diffy[svec], na.rm = T)
        }
        
        
        
        #historical
        festive_start <- (which(festive) - 1)*period
        
        for(start in festive_start) {
          festive_ind <- c(festive_ind, (1:period) + festive_start)
        }
        
        cell_F_historical_mean[[node]][kpi] <- mean(data[festive_ind], na.rm = T)
        cell_F_historical_sd[[node]][kpi] <- sd(data[festive_ind], na.rm = T)
        cell_historical_mean[[node]][kpi] <- mean(data[-festive_ind], na.rm = T)
        cell_historical_sd[[node]][kpi] <- sd(data[-festive_ind], na.rm = T)
        
        times <- ts(data[-festive_ind], frequency = period)
        Ftimes <- ts(data[festive_ind], frequency = period)
        
      }
      else if(!is.null(preCalcs)) {
        
        cell_F_seasonality <- preCalcs$F_seasonality
        cell_F_seasonality_sd <- preCalcs$F_season_sd
        cell_seasonality <- preCalcs$seasonality
        cell_seasonality_sd <- preCalcs$season_sd
        cell_F_historical_mean <- preCalcs$F_hist_mean
        cell_F_historical_sd$F_hist_sd
        cell_historical_mean$hist_mean
        cell_historical_sd$hist_sd       
        cell_diff_periodical_mean <- preCalcs$diff_season_mean
        cell_diff_periodical_sd <- preCalcs$diff_season_sd
        
      }
      
      for(hour in (0:(period-1))) {
        svec <- seq(1 + hour, length(data), period)
        festive_ind <- which(festive == T)
        #deseason data (residual)
        des_data[svec[-festive_ind]] <- des_data[svec[-festive_ind]] - cell_seasonality[[node]][kpi,hour + 1] 
        des_data[svec[festive_ind]] <- des_data[svec[festive_ind]] - cell_F_seasonality[[node]][kpi,hour + 1] 
      }
      
      if(calc_cdf == TRUE) {
        
        #calculate empirical kernel cdf from residual data
        if(length(unique(des_data)) < 2) {
          #Season is constant
          cell_cdf[[node]][[kpi]] <- list(x = c(des_data[1]),
                                          Fhat = c(1))
        }
        else {
          calced_cdf <- kCDF(des_data, na.rm = TRUE)
          cell_cdf[[node]][[kpi]] <- list(x = calced_cdf$x,
                                          Fhat = calced_cdf$Fhat)
        }
        
      }
      else if(!is.null(preCalcs)) {
        
        cell_cdf <- preCalcs$kernel_cdf
        
      }
      
      
    }
    
    
  }
  
  
  
  if(fits == TRUE) {
    
    print("Computing method fits...")  
    
    for(feature in 1:length(cellData)) {
      
      #what can't be regressed put to 0 - we probably don't care anyway
      cellData[[feature]][is.na(cellData[[feature]])] <- 0
      cellData[[feature]][is.nan(cellData[[feature]])] <- 0
      
    }  
    
    
    
    for(time in (1:final_timeslot)) {
      
      #check time slot
      curr_period <- ((time-1) %% period) + 1
      
      #check if festive
      curr_day <- floor((time-1)/period) + 1
      festivity <- festive[curr_day]
      
      print(paste("Time ",time))
      
      for(node in (1:node_number)) {
        curr_values <- vector("numeric", length = kpi_number)
        
        method_fit <- rep(list(vector("numeric", length = kpi_number)) ,method_number)
        
        
        ########################Parameter preparation
        
        #for each feature in the node, calculate various measures of unlikelihood
        for(kpi in (1:kpi_number)) {
          
          #print(paste("KPI", kpi))
          curr_values[kpi] <- cellData[[kpi]][node, time]
          
          if(festivity) {
            curr_seasonality <- cell_F_seasonality[[node]][kpi,]
            
          }
          else {
            curr_seasonality <- cell_seasonality[[node]][kpi,]
            
          }
          
          
          #find standardized values for distance and trend
          season_avg <- mean(curr_seasonality)
          season_rng <- (range(curr_seasonality)[2] 
                         - range(curr_seasonality)[1]) / 2
          
          if(season_rng < min_difference ) {
            season_rng <- min_difference
          }
          
          curr_val_std <- (curr_values[kpi] - season_avg) / season_rng 
          curr_season_std <- (curr_seasonality[curr_period] - season_avg) / season_rng 
          
          dist <- curr_val_std - curr_season_std
          
          if(curr_period > 1) {
            last_period <- curr_period - 1 
          }
          else {
            last_period <- period
          }          
          
          
          if(time >= 2) {
            
            
            prev_val_std <- (cell_memories[[node]][kpi,ncol(cell_memories[[node]])] - season_avg) / season_rng 
            prev_season_std <- (curr_seasonality[last_period] - season_avg) / season_rng 
            
            residual <- (curr_val_std - prev_val_std) - (curr_season_std - prev_season_std)
            
            
            last_descent_mean <- cell_descent_historical_mean[[node]][kpi]
            
            cell_descent_historical_mean[[node]][kpi] <- (cell_descent_historical_mean[[node]][kpi] * (time - 1) +
                                                            residual) / time
            cell_descent_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_descent_historical_sd[[node]][kpi] ^2) ) 
                                                             + (residual - last_descent_mean) * (residual - cell_descent_historical_mean[[node]][kpi])) / time)
            
            
            if(is.nan(cell_descent_historical_sd[[node]][kpi])) {
              cell_descent_historical_sd[[node]][kpi] <- 0
            }
            
            
          }
          
          if(time >= min_diff_window)  {
            
            last <- ncol(cell_memories[[node]])
            first <- last - min_diff_window + 2
            diff_memories <- c(cell_memories[[node]][kpi, first:last], curr_values[kpi])
            
            last <- curr_period
            first <- last - min_diff_window + 1
            if(first < 1) {
              
              window2 <- curr_seasonality[1:last]
              #check if previous window is festive or not
              first <- first + 24
              if(festive[curr_day -1]) {
                window1 <- cell_F_seasonality[[node]][kpi, first:period]
              }
              else {
                window1 <- cell_seasonality[[node]][kpi, first:period]
              }
              
              diff_season <- c(window1, window2)
            }
            else {
              diff_season <- curr_seasonality[first:last]
            }
            
            diff_memories <- (diff_memories - season_avg) / season_rng
            diff_season <- (diff_season - season_avg) / season_rng
            
            
            
            
            
            
            diff_residual <- (diff_memories - diff_season)
            
            myDiff <- diff(diff_residual, 
                           lag = diff_lag, 
                           differences = diff_differences)
            diff_diff <- myDiff[length(myDiff)]
            
            
            last_diff_mean <- cell_diff_historical_mean[[node]][kpi]
            
            cell_diff_historical_mean[[node]][kpi] <- (cell_diff_historical_mean[[node]][kpi] * (time - 1) +
                                                         diff_diff) / time
            cell_diff_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_diff_historical_sd[[node]][kpi] ^2) ) 
                                                          + (diff_diff - last_diff_mean) * (diff_diff - cell_diff_historical_mean[[node]][kpi])) / time)
            
            #             
            #             last_diff_mean_s <- cell_diff_periodical_mean[[node]][kpi, curr_period]
            #             
            #             cell_diff_periodical_mean[[node]][kpi, curr_period] <- (cell_diff_periodical_mean[[node]][kpi, curr_period] * (curr_day - 1) +
            #                                                          diff_diff) / curr_day
            #             cell_diff_periodical_sd[[node]][kpi, curr_period] <- sqrt((((curr_day - 1) * (cell_diff_periodical_sd[[node]][kpi] ^2) ) 
            #                                                           + (diff_diff - last_diff_mean_s) * (diff_diff - cell_diff_periodical_mean[[node]][kpi, curr_period])) / curr_day)
            #             
            #             
          }
          
          last_distance_mean <- cell_distance_historical_mean[[node]][kpi]
          
          cell_distance_historical_mean[[node]][kpi] <- (cell_distance_historical_mean[[node]][kpi] * (time - 1) +
                                                           dist) / time
          cell_distance_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_distance_historical_sd[[node]][kpi] ^2) ) 
                                                            + (dist - last_distance_mean) * (dist - cell_distance_historical_mean[[node]][kpi])) / time)
          
          
          if(is.nan(cell_distance_historical_sd[[node]][kpi])) {
            cell_distance_historical_sd[[node]][kpi] <- 0
          }
          
          
          if(cell_historical_sd[[node]][kpi] < min_difference) {
            cell_historical_sd[[node]][kpi] <- min_difference
          }
          
          
          ######################### Method application
          
          for (meth in 1:method_number) {
            
            if(time >= startMethod[meth] && selectMethod[meth]) {
              
              if(nameMethod[meth] == "season") {
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureSeasonFitness(curr_seasonality,
                                                                curr_values[kpi],
                                                                curr_period, replicas = 3)
                
              }            
              else if(nameMethod[meth] == "dist" ) {
                
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureDistanceThresholdPVal(curr_season_std,
                                                                        curr_val_std,
                                                                        cell_distance_historical_sd[[node]][kpi],
                                                                        threshold_confidence)
                
              }
              else if(nameMethod[meth] == "trend") {
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureTrendThresholdPVal(curr_season_std,
                                                                     prev_season_std,
                                                                     curr_val_std,
                                                                     prev_val_std,
                                                                     cell_descent_historical_sd[[node]][kpi],
                                                                     threshold_confidence) 
              }
              else if(nameMethod[meth] == "last") {
                
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureLastValuesFitness(cell_memories[[node]][kpi,],
                                                                    curr_values[kpi])   
              }
              else if(nameMethod[meth] == "history") {
                
                #print(nameMethod[meth])
                if(festivity) {
                  method_fit[[meth]][kpi] <- featureHistoryFitness(cell_F_historical_mean[[node]][kpi],
                                                                   cell_F_historical_sd[[node]][kpi],
                                                                   curr_values[kpi])
                  
                }
                else {
                  #print(nameMethod[meth])
                  method_fit[[meth]][kpi] <- featureHistoryFitness(cell_historical_mean[[node]][kpi],
                                                                   cell_historical_sd[[node]][kpi],
                                                                   curr_values[kpi])
                  
                }
                
              }
              else if(nameMethod[meth] == "cdf") {
                
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureResidualCDFPVal(curr_seasonality[curr_period] - curr_values[kpi],
                                                                  cell_cdf[[node]][[kpi]])
                
                
              }
              else if(nameMethod[meth] == "diff") {
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureDiffPVal(diff_season,
                                                           diff_memories,
                                                           cell_diff_historical_sd[[node]][kpi],
                                                           diff_thr_conf,
                                                           diff_lag,
                                                           diff_differences)
                
              }
              else if(nameMethod[meth] == "diff_season") {
                #print(nameMethod[meth])
                method_fit[[meth]][kpi] <- featureDiffPVal(diff_season,
                                                           diff_memories,
                                                           cell_diff_periodical_sd[[node]][kpi,curr_period],
                                                           diff_thr_conf,
                                                           diff_lag,
                                                           diff_differences)
                
              }
              
            }
            
            
            if(!is.null(method_fit[[meth]][[kpi]])){
              methFits[[node]][[kpi]][[meth]] <- c(methFits[[node]][[kpi]][[meth]], 
                                                   method_fit[[meth]][[kpi]])
            }
            
          }
          
          
        }
        
        #remember values
        cell_memories[[node]] <- cbind(cell_memories[[node]][,2:ncol(cell_memories[[node]])], curr_values)
        
        
      }
      
      
      
      ## END NODE OPERATION ##
      
    }
    
  }
  else if(!is.null(preCalcs)) {
    
    methFits <- preCalcs$method_fits
    cell_descent_historical_mean <- preCalcs$descent_mean
    cell_descent_historical_sd <- preCalcs$descent_sd
    cell_diff_historical_mean <- preCalcs$diff_mean
    cell_diff_historical_sd <- preCalcs$diff_sd
    #     cell_diff_periodical_mean <- preCalcs$diff_season_mean
    #     cell_diff_periodical_sd <- preCalcs$diff_season_sd
    cell_distance_historical_mean <- preCalcs$distance_mean
    cell_distance_historical_sd <- preCalcs$distance_sd
    
    
  }
  
  ######### Method weighting
  
  if(weighting == TRUE){
    print("Computing weights...")
    method_weight <- meth_weighting(methFits)
  }
  else if(!is.null(preCalcs)){
    method_weight  <- preCalcs$method_weights
    
  }
  else {
    method_weight <- 0
  }
  print("end")
  
  tr_result <- list(seasonality = cell_seasonality,
                    F_seasonality = cell_F_seasonality,
                    season_sd = cell_seasonality_sd,
                    F_season_sd = cell_F_seasonality_sd,
                    hist_mean = cell_historical_mean,
                    hist_sd = cell_historical_sd,
                    F_hist_mean = cell_F_historical_mean,
                    F_hist_sd = cell_F_historical_sd,
                    descent_mean = cell_descent_historical_mean,
                    descent_sd = cell_descent_historical_sd,
                    diff_mean = cell_diff_historical_mean,
                    diff_sd = cell_diff_historical_sd,
                    diff_season_mean = cell_diff_periodical_mean,
                    diff_season_sd = cell_diff_periodical_sd,
                    distance_mean = cell_distance_historical_mean,
                    distance_sd = cell_distance_historical_sd,
                    kernel_cdf = cell_cdf,
                    method_weights = method_weight,
                    method_fits = methFits)
  
  return(tr_result)
  
}



############################# Priority calculation methods for ensemble ##############


#returns: the pvalue of the value given the seasonality
#replicas: the function replicates the seasonality for a number of time to achieve better regression model (danger of overfitting!)
#max replicas : 5
featureSeasonFitness <- function(feat_cell_seasonality, feat_curr_value, curr_period, replicas = 2, fit = 6, visual = FALSE) {
  
  maxReplicas <- 5
  
  #shift cell seasonality to have last slot the previous value
  
  #feat_cell_seasonality = [1,2,3, ... ,23,24]
  #curr_period = 8
  # ---> previous_seasonality <- [8,9,10, ... ,6,7]
  
  previous_seasonality <- shiftLeft(feat_cell_seasonality, num = (curr_period - 1))
  
  library(forecast)
  
  replicas <- min(replicas, maxReplicas)
  
  if(abs(max(previous_seasonality) - min(previous_seasonality)) < 0.0001) {
    #nothing to predict, all equal
    
    fcasted <- max(previous_seasonality)
    sdev <- 0.0001
    
  }
  else {
    #fit an arima to the historical data
    #aa <- auto.arima(rep(previous_seasonality, replicas))
    aa <- Arima(rep(previous_seasonality, replicas), c(3,0,0), c(0,1,1), method = "CSS")
    
    ff <- forecast(aa, h = 1, level = c(68.27, 95.45))
    
    #find z-score (number of standard deviation of distance from forecast value)
    fcasted <- ff$mean[1]
    sdev <- abs(ff$lower[1] - ff$lower[2])
    
  }
  zscore <- (feat_curr_value - fcasted) / sdev
  
  #  pval <- 2*pnorm(-abs(zscore))
  pval <- pnorm(zscore)
  
  # (1 - pval) = minimum prediction interval threshold for given value, value within 1-pval most probable
  # pval = value is predicted within the pval less probable outcom
  
  if(visual == TRUE) {
    
    old.par <- par(mfrow = c(2,1), mai = c(0.6,1,0.5,0.5))
    
    plot(c(previous_seasonality, NA), type ="o", col ="blue", main = "Seasonal curve", pch = 21)
    points(length(previous_seasonality) + 1, feat_curr_value, col = "red")
    
    plot(ff, main = "Forecast")
    points(length(previous_seasonality) * replicas + 1, feat_curr_value, col = "red")
    
    par(old.par)
    
  }
  
  return (pval)  
  
}

#returns: the pvalue of the value given the seasonality
featureLastValuesFitness <- function(feat_cell_past_values, feat_curr_value, visual = FALSE, 
                                     order = c(0,1,0), season_order = c(0,1,0), prevModel = NULL) {
  
  #assumes last value in feat_cell_past_values is at the end
  
  library(forecast)
  
  if(abs(max(feat_cell_past_values) - min(feat_cell_past_values)) < 0.0001) {
    #nothing to predict, all equal
    
    fcasted <- max(feat_cell_past_values)
    sdev <- 0.0001
    
  }
  else {
    
    
    #fit an arima to the historical data
    
    if(is.null(prevModel)) {
    aa <- Arima(feat_cell_past_values, order, season_order, method = "CSS-ML")
    }
    else {
    }
      
    ff <- forecast(aa, h = 1, level = c(68.27, 95.45))
    
    #find z-score (number of standard deviation of distance from forecast value)
    fcasted <- ff$mean[1]
    sdev <- abs(ff$lower[1] - ff$lower[2])
    
  }
  
  zscore <- (feat_curr_value - fcasted) / sdev
  
  #  pval <- 2*pnorm(-abs(zscore))
  pval <- pnorm(zscore)
  
  
  # (1 - pval) = minimum prediction interval threshold for given value, value within 1-pval most probable
  # pval = value is predicted within the pval less probable outcome
  
  if(visual == TRUE) {
    
    old.par <- par(mfrow = c(2,1), mai = c(0.6,1,0.5,0.5))
    
    plot(c(feat_cell_past_values, NA), type ="o", col ="blue", main = "Past curve")
    points(length(feat_cell_past_values) + 1, feat_curr_value, col = "red")
    
    plot(ff, main = "Forecast")
    points(length(feat_cell_past_values) + 1, feat_curr_value, col = "red")
    
    par(old.par)
    
  }
  
  return (pval)  
  
}


#returns: the pvalue of the value given the historical mean and sd for the period
featureHistoryFitness <- function(feat_historical_period_mean, feat_historical_period_sd, feat_curr_value, visual = FALSE) {
  
  
  #find z-score (number of standard deviation of distance from historical value)
  zscore <- (feat_curr_value - feat_historical_period_mean) / feat_historical_period_sd
  
  #  pval <- 2*pnorm(-abs(zscore))
  pval <- pnorm(zscore)
  
  
  # (1 - pval) = minimum prediction interval threshold for given value, value within 1-pval most probable
  # pval = value is predicted within the pval less probable outcom
  
  if(visual == TRUE) {
    
    yminlim <- min(feat_historical_period_mean - (3 * feat_historical_period_sd),
                   feat_curr_value) - 0.5
    ymaxlim <- max(feat_historical_period_mean + (3 * feat_historical_period_sd),
                   feat_curr_value) + 0.5
    
    plot(1, type ="n", main = "Historical data", ylim = c(yminlim, ymaxlim))
    points(feat_curr_value, col = "black", pch = 16)
    
    abline(h = feat_historical_period_mean, col = "green")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean + 0.15, label = "Historical Mean", col = "green")
    
    abline(h = feat_historical_period_mean + feat_historical_period_sd, col = "blue")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean + feat_historical_period_sd + 0.15, label = "+1 sd", col = "blue")
    abline(h = feat_historical_period_mean - feat_historical_period_sd, col = "blue")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean - feat_historical_period_sd + 0.15, label = "-1 sd", col = "blue")
    
    abline(h = feat_historical_period_mean + 2 * feat_historical_period_sd, col = "purple")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean + 2 * feat_historical_period_sd + 0.15, label = "+2 sd", col = "purple")
    abline(h = feat_historical_period_mean - 2 *feat_historical_period_sd, col = "purple")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean - 2 * feat_historical_period_sd + 0.15, label = "-2 sd", col = "purple")
    
    abline(h = feat_historical_period_mean + 3 * feat_historical_period_sd, col = "red")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean + 3 * feat_historical_period_sd + 0.15, label = "+3 sd", col = "red")
    abline(h = feat_historical_period_mean - 3 *feat_historical_period_sd, col = "red")
    text(adj = c(0, NA), x = 0.6, y = feat_historical_period_mean - 3 * feat_historical_period_sd + 0.15, label = "-3 sd", col = "red")
    
  }
  
  return (pval)  
  
}


featureDistanceThresholdPVal <- function(feat_cell_season_val, 
                                         feat_curr_value, 
                                         distance_threshold, 
                                         threshold_confidence, 
                                         visual = FALSE,
                                         dir = -1) {
  
  #threshold_confidence indicates it is expected that value falls
  #within distance_threshold with that confidence
  
  residual <- feat_curr_value - feat_cell_season_val
  
  #given confidence of 95.45%, distance_treshold is 2 sd
  #given confidence of x, distance_treshold is zscore(1 - x) * sd?
  
  #  sDev <- distance_threshold / qnorm(1 - (1 - threshold_confidence)/2 )
  sDev <- distance_threshold / qnorm(1 - (1 - threshold_confidence))
  
  if(sDev == 0) {
    if(residual > 0) {
      zscore <- 5
    }
    else if(residual == 0) {
      zscore <- 0
    }
    else if(residual < 0) {
      zscore <- -5
    }
  }
  else{
    zscore <- residual / sDev
  }
  
  if(dir < 0) {
    #we care if zscore is low
    pval <- pnorm(zscore)
    
  }
  else if(dir > 0) {
    #we care if zscore is high
    pval <- pnorm(-zscore)
    
  }
  else {
    #we care if zscore is low or high
    pval <- pnorm(-abs(zscore))*2
  }
  
  
  if(visual == TRUE) {
    
    
    yminlim <- min(feat_cell_season_val - distance_threshold,
                   feat_curr_value) - 0.5
    ymaxlim <- max(feat_cell_season_val + distance_threshold,
                   feat_curr_value) + 0.5
    
    plot(1, type ="n", main = "Distance threshold", ylim = c(yminlim, ymaxlim))
    rect(par("usr")[1], feat_cell_season_val - distance_threshold, 
         par("usr")[2], feat_cell_season_val + distance_threshold, col = "whitesmoke")
    points(feat_cell_season_val, col = "blue", pch = 21)
    points(feat_curr_value, col = "black", pch = 16)
    
    abline(h = feat_cell_season_val + distance_threshold, col = "red")
    abline(h = feat_cell_season_val - distance_threshold, col = "red")
    
  }
  
  
  return(pval)
  
}

featureDiffPVal <- function(season_window, measure_window,
                            threshold, threshold_confidence, 
                            lag, difference,
                            dir = -1) {
  
  #threshold_confidence indicates it is expected that value falls
  #within distance_threshold with that confidence
  
  season_residual <- measure_window - season_window
  
  myDiff <- diff(season_residual, lag = lag, differences = difference)
  residual <- myDiff[length(myDiff)]
  
  #given confidence of 95.45%, distance_treshold is 2 sd
  #given confidence of x, distance_treshold is zscore(1 - x) * sd?
  
  # sDev <- descent_threshold / qnorm(1 - (1 - threshold_confidence)/2 )
  sDev <- threshold / qnorm(threshold_confidence)  
  
  if(sDev == 0) {
    if(residual > 0) {
      zscore <- 5
    }
    else if(residual == 0) {
      zscore <- 0
    }
    else if(residual < 0) {
      zscore <- -5
    }
  }
  else{
    zscore <- residual / sDev
  }
  
  if(dir < 0) {
    #we care if zscore is low
    pval <- pnorm(zscore)
    
  }
  else if(dir > 0) {
    #we care if zscore is high
    pval <- pnorm(-zscore)
    
  }
  else {
    #we care if zscore is low or high
    pval <- pnorm(-abs(zscore))*2
  }
  
  
  return(pval)
  
  
}


featureTrendThresholdPVal <- function(feat_cell_season_val, feat_prev_cell_season_val,
                                      feat_curr_value, feat_prev_value, 
                                      descent_threshold, 
                                      threshold_confidence,
                                      visual = FALSE,
                                      dir = -1) {
  
  #threshold_confidence indicates it is expected that value falls
  #within distance_threshold with that confidence
  
  curr_residual <- feat_curr_value - feat_prev_value
  season_residual <- feat_cell_season_val - feat_prev_cell_season_val
  
  residual <- curr_residual - season_residual
  
  #given confidence of 95.45%, distance_treshold is 2 sd
  #given confidence of x, distance_treshold is zscore(1 - x) * sd?
  
  # sDev <- descent_threshold / qnorm(1 - (1 - threshold_confidence)/2 )
  sDev <- descent_threshold / qnorm(threshold_confidence)  
  
  if(sDev == 0) {
    if(residual > 0) {
      zscore <- 5
    }
    else if(residual == 0) {
      zscore <- 0
    }
    else if(residual < 0) {
      zscore <- -5
    }
  }
  else{
    zscore <- residual / sDev
  }
  
  if(dir < 0) {
    #we care if zscore is low
    pval <- pnorm(zscore)
    
  }
  else if(dir > 0) {
    #we care if zscore is high
    pval <- pnorm(-zscore)
    
  }
  else {
    #we care if zscore is low or high
    pval <- pnorm(-abs(zscore))*2
  }
  
  if(visual == TRUE) {
    
    old.par <- par(mfrow = c(2,1), mai = c(0.6,1,0.5,0.5))
    
    yminlim <- min(feat_cell_season_val, feat_prev_cell_season_val,
                   feat_curr_value, feat_prev_value) - 0.5
    ymaxlim <- max(feat_cell_season_val, feat_prev_cell_season_val,
                   feat_curr_value, feat_prev_value) + 0.5
    
    plot(1, type ="n", main = "Trends", xlim = c(-0.5, +1.5), ylim = c(yminlim, ymaxlim))
    points(c(0, 1), c(feat_prev_cell_season_val, feat_cell_season_val), col = "blue", pch = 21)
    lines(c(0, 1), c(feat_prev_cell_season_val, feat_cell_season_val), col = "blue")
    points(c(0, 1), c(feat_prev_value, feat_curr_value), col = "black", pch = 16)
    lines(c(0, 1), c(feat_prev_value, feat_curr_value), col = "black")
    
    plot(1, type ="n", main ="Threshold", xlim = c(-1, 1), ylim = c(-1, 1))
    abline(a = 0, b = (feat_cell_season_val - feat_prev_cell_season_val), col = "blue")
    abline(a = 0, b = (feat_curr_value - feat_prev_value), col = "black")
    abline(a = 0, b = (feat_cell_season_val - feat_prev_cell_season_val) + descent_threshold, col = "red")
    abline(a = 0, b = (feat_cell_season_val - feat_prev_cell_season_val) - descent_threshold, col = "red")   
    #further sd thresholds
    abline(a = 0, b = (feat_cell_season_val - feat_prev_cell_season_val) + 2 * descent_threshold, col = "purple")
    abline(a = 0, b = (feat_cell_season_val - feat_prev_cell_season_val) - 2 * descent_threshold, col = "purple")    
    
    
    par(old.par)
    
  }
  
  
  return(pval)
  
}

#Generic distribution replacement for distance from the mean
#Right now Only uni-directional decrease
featureResidualCDFPVal <- function(curr_value, cdf) {
  
  
  bef <- max(which(cdf$x <= curr_value))
  aft <- min(which(cdf$x >= curr_value))
  
  if(bef < 1) {
    p_val <- 0
  }
  else if (aft > length(cdf$x)) {
    p_val <- 1
  }
  else {
    
    p_dist <- cdf$Fhat[aft] - cdf$Fhat[bef]
    x_dist <- cdf$x[aft] - cdf$x[bef]
    
    if(is.nan(x_dist) || is.na(x_dist)) {
      
      print("debud")
    }
    
    if(x_dist != 0) {
      
      nearB <- (curr_value - cdf$x[bef]) / x_dist
      
      p_val <- cdf$Fhat[bef] + p_dist * (nearB)
      
    }
    else {
      
      p_val <- cdf$Fhat[bef]
    }
    
  }
  
  return(p_val)
  
}


############################# Fisher Test and priority combination ##############


fastGeneralFisherScan <- function(priorities, anomThreshold = 0.95) {
  
  
  #order priorities
  prior_order <- sort(priorities, decreasing = TRUE,index.return = T)$ix
  
  subset <- NULL
  Fun <- rep( 0  , length(prior_order))
  
  #in order, calculate F(S)
  for (size in c(1:length(prior_order))) {
    
    nodes <- prior_order[1:size]
    
    Fun[size] <- fisherCombinedPriority(1 - priorities[nodes])
    
    #Fun[size:length(Fun)] <- Fun[size:length(Fun)] - (2 * log(1 - priorities[node]))
    
  }
  
  string <- NULL
  
  final_scores <- Fun
  
  #Best result
  
  bestSize <- 0
  bestScore <- 0
  
  for(size in c(1:length(prior_order))) {
    if(final_scores[size] >= bestScore) {
      bestScore <- final_scores[size]
      bestSize <- size  
    }
  }
  for(size in c(1:bestSize)) {
    string <- paste(string," ", prior_order[size])
  }
  
  if(bestScore > anomThreshold) {
    #print(paste("Set ",size,", neighbours (",string,"), score ",bestScore))
  }
  
  return(list(vote = bestScore,
              nodes = prior_order[1:bestSize]))
  
}

fastGeneralWeightedFisherScan <- function(priorities, weights, anomThreshold = 0.95, base_priority = 0.68,
                                          best = TRUE) {
  
  
  #baselines <- fisher.baselineDiscovery(priorities, weights, baseStart = base_priority, iterations = 10)
  baselines <- rep(0.6217, length(priorities))
  
  #transform priorities according to their weight
  for(p in (1:length(priorities))) {
    
    priorities[p] <- priorityWeight(priorities[p], weights[p], baselines[p])
    
  }
  
  #order priorities
  prior_order <- sort(priorities, decreasing = TRUE,index.return = T)$ix
  
  subset <- NULL
  Fun <- rep( 0  , length(prior_order))
  
  #in order, calculate F(S)
  for (size in c(1:length(prior_order))) {
    
    nodes <- prior_order[1:size]
    
    Fun[size] <- fisherCombinedPriority(1 - priorities[nodes])
    
    #Fun[size:length(Fun)] <- Fun[size:length(Fun)] - (2 * log(1 - priorities[node]))
    
  }
  
  string <- NULL
  
  final_scores <- Fun
  
  if(best) {
    #Best result
    
    bestSize <- 0
    bestScore <- 0
    
    
    for(size in c(1:length(prior_order))) {
      if(final_scores[size] >= bestScore) {
        bestScore <- final_scores[size]
        bestSize <- size  
      }
    }
    
    for(size in c(1:bestSize)) {
      string <- paste(string," ", prior_order[size])
    }
    
    
    return(list(vote = bestScore,
                nodes = prior_order[1:bestSize]))
  }
  else {
    #Biggest over threshold
    
    bestScore <- anomThreshold
    bestSize <- 0
    
    for(size in c(1:length(prior_order))) {
      if(final_scores[size] >= anomThreshold) {
        bestScore <- final_scores[size]
        bestSize <- size
      }
    }
    
    
    return(list(vote = bestScore,
                nodes = prior_order[1:bestSize]))
    
  }
}


#pvector: pvalues for the univariate\multivariate case
#returns: the priority of the combined
priority_function <- function(pvector, alfa) {
  
  
  #return the fisher combined probability -> the priority function for itself
  #this is thanks to its chi-square distribution, the sum of two sample of which results in
  #a chi-square distribution
  
  #this property allows fisher's test to satisfy the strong LTSS property
  #it is preferable to common ltss count-based function that do not account for
  #multiple hypotesis testing
  return(fisherCombinedPriority(pvector))
  
  
}

fisher.show <- function(priorities) {
  
  changingp <- c(1000:0)/1000
  res <- NULL
  for(i in 1:length(changingp)) {
    res <- c(res, fisherCombinedPriority(c(1 - priorities, changingp[i])))
  }
  sting <- paste(priorities, collapse = ", ")
  string <- paste("Fisher curve, starting priorities", sting)
  #plot(1 - changingp, 1 - changingp, type = "l", main = string)
  plot(1 - changingp, res, col = "red", type = "l", main = string, ylim = c(0,1))
  abline(h = fisherCombinedPriority(c(1 - priorities)), col = "blue")
  
  lineL <- 10
  for(j in 1:length(priorities)) {
    points((1:lineL)/1000, rep(priorities[j],lineL), col = "purple", type = "l")
  }
  
}
fisher.show2 <- function(normPriorities, changPriority, changWeight) {
  
  changingb <- c(1000:0)/1000
  res <- NULL
  for(i in 1:length(changingb)) {
    res <- c(res, fisherCombinedPriority(c(1 - normPriorities, 
                                           1 - priorityWeight(changPriority, changWeight, changingb[i]))))
  }
  
  plot(res, type = "o")
}

#normPriorities already weighted with their baselines
fisher.baseChange <- function(normPriorities, minb = 0, maxb = 1, granularity = 100) {
  
  changingb <- c((minb*granularity):(maxb*granularity))/granularity
  res <- NULL
  minres <- Inf
  minb <- -1
  for(i in 1:length(changingb)) {
    
    orig <- fisherCombinedPriority(c(1 - normPriorities))
    
    add <- fisherCombinedPriority(c(1 - normPriorities,  1 - changingb[i]))
    
    if(abs(add - orig) < minres) {
      minres <- abs(add-orig)
      minb <- changingb[i]
    }
    else if (abs(add - orig) == minres) {
      
      minb <- c(minb, changingb[i])
      
    }
    res <- c(res, add - orig)
  }
  #plot(abs(res), type = "l")
  
  return(list(base = minb,
              error = minres))
}


fisher.baselineTest <- function(restarts = 100) {
  
  baselines <- NULL
  
  for(i in 1:restarts) {
    
  print(paste("Restart ", i))
  priorities <- runif(30, 0, 1)
  weights <- runif(30, 0, 1)
  baseStart = runif(1, 0, 1)
  
  minb = 0.53
  maxb = 0.73
  granularity = 1000
  
  
  baselines <- rbind(baselines, fisher.baselineDiscovery(priorities, weights, baseStart, 
                                                         iterations = 10,
                                                         minb, maxb, granularity))
  
  }
  
  return(baselines)
}


fisher.baselineErrorTest <- function(restarts = 100) {
  
  errorsmin <- NULL
  errorsmax <- NULL
  
  minb = 0.584
  maxb = 0.684
  
  for(i in 1:restarts) {
    
    print(paste("Restart ", i))
    priorities <- runif(30, 0, 1)
    weights <- runif(30, 0, 1)
    
    weightPriorities <- vector("numeric", length(priorities))
    
    #min
    baselines <- rep(minb, length(priorities))
    for(i in 1:length(priorities)) {
      #compute weightede priority with old weights
      weightPriorities[i] <- priorityWeight(priorities[i], weights[i], baselines[i])
      
    }
    for(i in  1:length(priorities)) {
      
      orig <- fisherCombinedPriority(c(1 - weightPriorities[-i]))
      
      add <- fisherCombinedPriority(c(1 - weightPriorities[-i],  1 - baselines[i]))
      
      errorsmin <- c(errorsmin, abs(add - orig))
    }
    
    #max
    baselines <- rep(maxb, length(priorities))
    for(i in 1:length(priorities)) {
      #compute weightede priority with old weights
      weightPriorities[i] <- priorityWeight(priorities[i], weights[i], baselines[i])
      
    }
    for(i in  1:length(priorities)) {
      
      orig <- fisherCombinedPriority(c(1 - weightPriorities[-i]))
      
      add <- fisherCombinedPriority(c(1 - weightPriorities[-i],  1 - baselines[i]))
      
      errorsmax <- c(errorsmax, abs(add - orig))
    }
    
      
  }
  
  return(abs(errorsmax - errorsmin))
}

fisher.baselineDiscovery <- function(priorities, weights, 
                                  baseStart = 0.6, iterations = 10, minb = 0,
                                  maxb = 1, granularity = 100) {
  
  
  baselines <- rep(baseStart, length(priorities))
  
  errorGrowth <- NULL
  histErrors <- NULL
#   
#   string <- paste(baselines, collapse = ", ")
#   print(paste("Start baselines: ",string))
#   
#   
  lasterrors <- NULL
  
  weightPriorities <- priorities
  for(i in 1:length(priorities)) {
    #compute weightede priority with old weights
    weightPriorities[i] <- priorityWeight(priorities[i], weights[i], baselines[i])
    
  }
  
  #Initial error
  for(i in  1:length(priorities)) {
    
    orig <- fisherCombinedPriority(c(1 - weightPriorities[-i]))
    
    add <- fisherCombinedPriority(c(1 - weightPriorities[-i],  1 - baselines[i]))
    
    lasterrors <- c(lasterrors, abs(add - orig))
  }
  
  
  #Iterative process
  for(it in 1:iterations) {

    newBaselines <- baselines
    weightPriorities <- priorities
    
    for(i in 1:length(priorities)) {
      #compute weightede priority with old weights
      weightPriorities[i] <- priorityWeight(priorities[i], weights[i], baselines[i])
      
    }
    
    errors <- NULL
    for(i in 1:length(priorities)) {
      #compute baselines that minimizes difference
      fisher_test <- fisher.baseChange(weightPriorities[-i], 
                                       minb, maxb, granularity)
      newBaselines[i] <- fisher_test$base 
      errors <- c(errors, fisher_test$error)
    }
    
    errorGrowth <- errors - lasterrors #c(errorGrowth, sum(errors ^ 2))
    
    newDiff <- baselines - newBaselines
    
    growing <- which(errorGrowth > 0)
    
    newBaselines[growing] <- baselines[growing] 
    
    
    baselines <- newBaselines
    lasterrors <- errors
  }
  
return(baselines)
  
}

#input: the pvalues
#Returns fisher combined test (higher = worse!) a priority
fisherCombinedPriority <- function(pvector) {
  
  max_priority <- 0.99999
  
  #to avoid problems, pvalues = 0 are considered with small but not zero value
  pvector[pvector < (1 - max_priority)] <- (1 - max_priority)
  
  #Fisher's combined probability test
  #
  fisher <- pchisq(-2 * (sum(log(pvector))), df = 2*length(pvector))
  
  return(fisher)
  
}

#weight between 0 and 1
#low weight will be indicated via base_priority
# everything with higher than base_priority will be lowered to the base_priority if their weight is low
# everything with lower than base_priority will be raised to the base_priority if their weight is low

#Base priority represents the priority of a node that doesn't give me any information on whether to 
#accept or reject the normality\anomalousness hypothesis
priorityWeight <- function(priority, weight, base_priority) {
  
  pval <- 1 - priority
  base_confidence <- 1 - base_priority
  
  myPower <- pwr( (pval - base_confidence), 1 / weight) + base_confidence
  dist <- myPower - pval
  
  if(pval < base_confidence) {
    wpval <- myPower - (dist * (abs(pval - base_confidence)/base_confidence)^(1/weight))
  }
  else {
    wpval <- myPower - (dist * (abs(pval - base_confidence)/ (1 - base_confidence) )^(1/weight))
    
  }
  
  return(1 - wpval)
  
}


#input: weights, and the pvalues
#Returns fisher combined test (higher = worse!) a priority
fisherWeightedCombinedPriority <- function(weights, pvector) {
  
  if(max(weights) > 1) {
    weights <- weights / max(weights)
  }
  
  priorities <- 1 - pvector
  #wpriorities <- weights * priorities
  wpriorities <- pchisq(qgamma(priorities, shape = weights / 2, scale = 2), df = 1)
  wpvalues <- 1 - wpriorities
  
  #Fisher's combined probability test
  #
  fisher <- pchisq(-2 * (sum(log(wpvalues))), df = 2*length(pvector))
  
  return(fisher)
  
}


#input: weights and pvalues
#weights are normalized to sum to 1 before applying
#Returns lancaster generalized fisher test, a priority
lancasterGeneralizedFisher <- function(weights, pvalues) {
  
  weights <- weights / sum(weights)
  
  priorities <- 1 - pvalues
  
  lancaster <- pchisq( sum(qgamma(priorities, shape = weights / 2, scale = 2)),
                       df = sum(weights))
  
  return(lancaster)
  
}

testLancaster <- function() {
  
  lancaster_growth <- NULL
  other <- NULL
  other2 <- NULL
  other3 <- NULL
  other4 <- NULL
  for(i in 1:100) {
    
    #assume growing weights*priority
    i <- i/100
    
    pi <- runif(1,0,1)
    wi <- i / pi
    
    
    
    other <- c(other, qgamma(pi, shape = wi / 2, scale = 2))
    other2 <- c(other2, qgamma(pi, shape = 1 / 2, scale = 2))
    
    other3 <- c(other3, pchisq(sum(qgamma(pi, shape = wi / 2, scale = 2)),
                               df = sum(wi)))
    other4 <- c(other4, pchisq(sum(qgamma(pi, shape = 1 / 2, scale = 2)),
                               df = sum(1)))
    
    #thus we have wi*pi = i, therefore growing at each iteration
    #let's check that wi*pi can be a valid priority metric (combined priority raises with it - convexity)
    
    lancaster_growth <- c(lancaster_growth,
                          lancasterGeneralizedFisher(c(0.5, wi), c(0.7, pi)))
    
    
  }
  
  #various ordering
  ord <- sort(other, index.return = T)$ix
  ord2 <- sort(other2, index.return = T)$ix
  ord3 <- sort(other3, index.return = T)$ix
  ord4 <- sort(other4, index.return = T)$ix
  
  #plotting this shows that sets containing couples wi, pi with high priority do not always have higher
  #values of lanc.gf than sets containing lower couples.
  
  return(lancaster_growth)
  
}

testLancaster2 <- function() {
  
  lancaster_growth <- NULL
  gamma_growth <- NULL
  
  for(i in 1:100) {
    
    #generate random weight and pvalues
    
    wi <- runif(1,0,1)
    pi <- runif(1,0,1)
    
    #generate quality (G(s))
    gamma_growth <- c(gamma_growth, qgamma(1 - pi, shape = wi / 2, scale = 2))
    
    lancaster_growth <- c(lancaster_growth,
                          lancasterGeneralizedFisher(c(0.5, wi), c(0.7, pi)))
    
    
  }
  
  #let's check if ordering lancaster_growth according to gamma_growth yields results
  
  sortorder <- sort.int(gamma_growth, index.return = T)$ix
  lancaster_ord <- lancaster_growth[sortorder]
  
  return(lancaster_ord)
  
}

testWeightedFisher <- function() {
  
  fisher_growth <- NULL
  
  for(i in 1:100) {
    
    #assume growing weights*priority
    i <- i/100
    
    pi <- runif(1,0,1)
    wi <- i / pi
    
    #thus we have wi*pi = i, therefore growing at each iteration
    #let's check that wi*pi can be a valid priority metric (combined priority raises with it - convexity)
    
    fisher_growth <- c(fisher_growth,
                       fisherCombinedPriority(c(0.6*0.8, 1 - wi*pi)))
    
    
    
  }
  
  return(fisher_growth)
  
}

testWeightedFisher2 <- function() {
  
  fisher_growth <- NULL
  weighted_growth <- NULL
  
  for(i in 1:1000) {
    
    #Calculate random weights and pvalues
    wi <- runif(1,0,1)
    pi <- runif(1,0,1)
    
    #Calculate priority from weight and pval
    #for a confidence threshold of 0.35
    confid <- 0.35
    powered <-  pwr((pi - confid) , 1 / wi) + confid
    weighted_growth <- c(weighted_growth, 1 - powered)
    
    fisher_growth <- c(fisher_growth,
                       fisherCombinedPriority(c(0.8, powered)))
    
  }
  
  #let's check if ordering lancaster_growth according to gamma_growth yields results
  
  sortorder <- sort.int(weighted_growth, index.return = T)$ix
  fisher_ord <- fisher_growth[sortorder]
  
  
  return(fisher_ord)
  
}




##################### Anomaly clusters operation ###########

anomaly_cluster.search_connection <- function(cluster, new_anomaly) {
  
  if(length(intersect(cluster$nodes, new_anomaly)) > 0 ) {
    
    return(TRUE)
    
  }
  else {
    return(FALSE)
  }
  
  
}

anomaly_cluster.merge <- function(anomaly_cluster, new_anomaly, new_priority) {
  
  merged <- FALSE
  mergedWith <- NULL
  
  if(!(is.null(anomaly_cluster) || length(anomaly_cluster) < 1)) {
    
    for(clust in 1:length(anomaly_cluster)) {
      
      if(anomaly_cluster.search_connection(anomaly_cluster[[clust]], new_anomaly) == TRUE) {
        
        anomaly_cluster[[clust]]$nodes <- union(anomaly_cluster[[clust]]$nodes, 
                                                new_anomaly)
        
        #SCORRETTO: ma per ora ok
        anomaly_cluster[[clust]]$priority <- mean(anomaly_cluster[[clust]]$priority, 
                                                  new_priority)
        
        merged <- TRUE
        mergedWith <- c(mergedWith, clust)
      }
      
    }
  }
  
  if(merged == FALSE) {
    #add to cluster
    anomaly_cluster <- c(anomaly_cluster, list(list(nodes = new_anomaly,
                                                    priority = new_priority)))
    
    
    
  }
  else if(length(mergedWith) >= 2) {
    
    #merge together all merged cluster
    fin_clust <- mergedWith[1]
    for(mer_clust in mergedWith[2:length(mergedWith)]) {
      
      anomaly_cluster[[fin_clust]]$nodes <- union(anomaly_cluster[[fin_clust]]$nodes, 
                                                  anomaly_cluster[[mer_clust]]$nodes)
      
      anomaly_cluster[[fin_clust]]$priority <- mean(anomaly_cluster[[fin_clust]]$priority, 
                                                    anomaly_cluster[[mer_clust]]$priority)
      
    }
    
    anomaly_cluster <- anomaly_cluster[-mergedWith[2:length(mergedWith)]]
    
  }
  
  return(anomaly_cluster)
  
}


##################### Support Functions #############


meth_weighting <- function(method_fits) {
  
  anomaly_significance <- 0.9545
  
  node_number <- length(method_fits)
  kpi_number <- length(method_fits[[1]])
  method_number <- length(method_fits[[1]][[1]])
  
  
  method_weight <- NULL
  
  for(node in 1:length(method_fits)) {
    print(paste("node ",node))
    new_mat <- NULL
    for(kpi in 1:length(method_fits[[node]])) {
      new_column <- NULL
      for(meth in 1:length(method_fits[[node]][[kpi]])) {
        
        new_weight <- single_weighting(method_fits[[node]][[kpi]][[meth]],
                                       1,
                                       anomaly_significance)
        
        new_column <- rbind(new_column, new_weight)
        
        
      }
      new_mat <- cbind(new_mat, new_column)
    }
    method_weight <- c(method_weight, list(new_mat))
  }
  
  return(method_weight)
}

single_weighting <- function(fit_series, start_weight,
                             anomaly_significance,
                             adjust = TRUE,
                             visual = FALSE) {
  
  
  
  methWeightP <-  vector("numeric", length = 0)
  
  time_slots <- length(fit_series)
  method_weight <- start_weight
  
  
  normal_rate <- length(which(1 - fit_series < anomaly_significance))/length(fit_series)
  anomaly_rate <- length(which(1 - fit_series > anomaly_significance))/ length(fit_series)
  
  
  anom_to_normal <- (anomaly_rate + (1 - anomaly_significance)) / normal_rate
  
  weight_up <- 1- anomaly_significance
  weight_down <- (weight_up * (anomaly_rate)/(1 - anomaly_significance))
  
  
  
  change <- NULL
  ##### Estimate
  
  an_error <- (anomaly_rate - (1 - anomaly_significance))*2.5
  method_weight <- min(max(0, 1 - an_error),1)
  adjust <- F
  
  
  for(time in 1:time_slots) {
    
    
    curr_fit <- fit_series[time]
    
    #CHECK THIS
    baseline <- 0.6
    
    
    if(!is.na(curr_fit)) {
      
      
      weightedP <- priorityWeight(1 - curr_fit,
                                  method_weight,
                                  #start_weight,
                                  baseline)
      
      methWeightP <- c(methWeightP,
                       weightedP)
      
      
      if(adjust == TRUE) {
        if(weightedP > anomaly_significance) {
          #punish
          method_weight <- max(0,method_weight - (weight_down))# * (weightedP - anomaly_significance)  / (1 - anomaly_significance)
          change <- c(change, - weight_down)
        }
        else if(weightedP < anomaly_significance) {
          #reward
          method_weight <- min(1, method_weight + weight_up) #* (anomaly_significance - weightedP) / anomaly_significance
          change <- c(change, + weight_up)
        }
      }
      
    }
    
  }
  
  method_weight <- max(min(1, method_weight),0)
  
  if(visual == TRUE) {
    if(adjust == TRUE) {
      
      methWeightP <-  vector("numeric", length = 0)
      
      
      for(time in 1:time_slots) {
        
        
        curr_fit <- fit_series[time]
        
        #CHECK THIS
        baseline <- 0.6
        
        
        if(!is.na(curr_fit)) {
          
          weightedP <- priorityWeight(1 - curr_fit,
                                      method_weight,
                                      baseline)
          
          methWeightP <- c(methWeightP,
                           weightedP)
          
          
        }
        
      }
    }
    
    #### Validation
    
    dens <- density(methWeightP)
    plot(dens, main = paste("Distribution of weighted priority with w = ",method_weight))
    
    plot(ecdf(methWeightP), main = paste("CDF of weighted priority with w = ",method_weight), xlim = c(0,1), ylim = c(0,1))
    #True Negatives
    rect(0, anomaly_significance, anomaly_significance, 1, col = "#FF003322")
    #False Positives
    rect(anomaly_significance, 0, 1, anomaly_significance, col = "#FF003322")
    #True Positives
    rect(anomaly_significance, anomaly_significance, 1,1, col = "#11FF2235")
    
    #points(ecdf(methWeightP))
    
    abline(v = anomaly_significance, col = "red")
    abline(h = anomaly_significance, col = "red")
  }
  
  return(method_weight)
}

#Computes discounted similarity between windows to compute weights between neighbors
similarity_check <- function(prev_priority_own, prev_priority_neigh) {
  
  discount <- 0.6
  
  dissimilar <- 0
  max_dissimilarity <- 0
  
  for(mem in (1:length(prev_priority_own))) {
    dissimilar <- dissimilar + abs(prev_priority_own[mem] - prev_priority_neigh[mem]) *discount ^(mem - 1)
    max_dissimilarity <- max_dissimilarity + 1 *discount ^(mem - 1)
  }
  
  return(1 - (dissimilar / max_dissimilarity))
  
}

#shift cells of a vector one position to the right, with wrapping
#shiftRight([1,2,3,4,5], num = 2) --> [4,5,1,2,3]
shiftRight <- function(vec, num = 1) {
  
  vec <- c(vec[(length(vec) - num + 1):length(vec)],vec[1:(length(vec)-num)])
  return(vec)
  
}

#shift cells of a vector num position to the left, with wrapping
#shiftLeft([1,2,3,4,5], num = 2) --> [3,4,5,1,2]
#num must be < length(vec)

shiftLeft <- function(vec, num = 1) {
  
  if(is.null(vec) || length(vec) < 1) {
    return(NULL)
  }
  
  while(num >= length(vec)) {
    num <- num - length(vec)
  }
  
  if(num > 0) {
    vec <- c(vec[(1 + num) : length(vec)],vec[1:num])
  }
  
  return(vec)
  
}


##################### Other Functions - Preprocess, visualization, analysis ###########

#Used from console
#Prints the day and hour relative to a time slot.
#Use offset to indicate number of starting day if the database does not start from 1
#ex: for week 1 - offset = 1;
#week 2 - offset = 9
#week 3 - offset = 16
#week 4 - offset = 23
col.printDate <- function(col, offset = 1) {
  
  print(paste("Day ", floor((col - 1)/24) + offset, ", hour ",(col -1)%%24,":00", collapse = ""))
  
  
}



subFeatureTrain <- function(train_data, kpis) {
  
  node_number <- length(train_data$seasonality)
  
  for(node in 1:node_number) {
    
    train_data$seasonality[[node]] <- matrix(train_data$seasonality[[node]][kpis,], nrow = length(kpis))
    train_data$F_seasonality[[node]] <- matrix(train_data$F_seasonality[[node]][kpis,] , nrow = length(kpis))
    train_data$season_sd[[node]] <- matrix(train_data$season_sd[[node]][kpis,] , nrow = length(kpis))
    train_data$F_season_sd[[node]] <- matrix(train_data$F_season_sd[[node]][kpis,] , nrow = length(kpis))
    train_data$hist_mean[[node]] <- train_data$hist_mean[[node]][kpis]
    train_data$hist_sd[[node]] <- train_data$hist_sd[[node]][kpis] 
    train_data$F_hist_mean[[node]] <- train_data$F_hist_mean[[node]][kpis] 
    train_data$F_hist_sd[[node]] <-  train_data$F_hist_sd[[node]][kpis] 
    train_data$descent_mean[[node]] <- train_data$descent_mean[[node]][kpis] 
    train_data$descent_sd[[node]] <- train_data$descent_sd[[node]][kpis] 
    train_data$diff_mean[[node]] <- train_data$diff_mean[[node]][kpis] 
    train_data$diff_sd[[node]] <- train_data$diff_sd[[node]][kpis]     
    train_data$diff_season_mean[[node]] <- matrix(train_data$diff_season_mean[[node]][kpis,], nrow = length(kpis))
    train_data$diff_season_sd[[node]] <- matrix(train_data$diff_season_sd[[node]][kpis,] , nrow = length(kpis))
    train_data$distance_mean[[node]] <- train_data$distance_mean[[node]][kpis] 
    train_data$distance_sd[[node]] <- train_data$distance_sd[[node]][kpis]
    train_data$kernel_cdf[[node]] <- train_data$kernel_cdf[[node]][kpis]
    train_data$method_weights[[node]] <- matrix(train_data$method_weights[[node]][,kpis], ncol = length(kpis))
    train_data$method_fits[[node]] <- train_data$method_fits[[node]][kpis]
    
    
  }
  
  return(train_data)
} 

subFeatureTest <- function(test_set, kpis) {
  
  return(test_set[kpis])
  
  
}

#Returns a list of data from the whole dataset
#Assumes the dataset first two column for each feature are the coordinates
#Does not reduce number of locations
subData <- function(dataset, fromDay, toDay, features = 1:length(dataset)) {
  
  sub <- NULL
  cols <- NULL
  for(i in 1:length(fromDay)) {
    
    cols <- c(cols, 2 + c(((fromDay[i] - 1)*24 + 1) : (toDay[i]*24)))
  }
  
  for(featInd in features) {
    
    feat <- dataset[[featInd]]
    
    sub <- c(sub, list(as.matrix(feat[,cols])))
    
  }
  
  return(sub)
  
}

#Returns a vector with T if a day is festive and F if it is not.
#The days refer to the november dataset, with 1 November being a Friday
subDays <- function(fromDay, toDay) {
  
  fres <- NULL
  
  for(i in 1:length(fromDay)) {
    
    
    days <- c(fromDay[i]:toDay[i])
    
    res <- rep(F, length(days))
    
    #if 1 November, T
    res[which(days == 1)] <- T
    
    #if saturday (2) or sunday (3)
    res[which(days%%7 >= 2 &  days%%7 <= 3)] <- T
    
    fres <- c(fres, res)
  }
  
  return(fres)
  
}

#Gets a table obtained from a csv, cleans it if not already done (set flag!)
#If needed, aggregates data according to coordinates
#Removes rows without data, imputes na
#
#After the process return the data.frame with x,y, and time wide values
#
#To use with other features, time sets, merge by x y to find common coordinates.
fullTreatment <- function(table, clean = FALSE, aggregateCoord = TRUE) {
  
  source("preprocess.R")
  
  if(clean){
    cleaned1 <- toClean1(table)
  }
  else {
    cleaned1 <- table
  }
  
  if(aggregateCoord) {
    library(reshape)
    print("Melting")
    melted <- toMelt1(cleaned1)
    print("Casting")
    casted <- cleanToNumeric(cast(melted, x + y ~ variable, fun = mean, na.rm = T))
    
  }
  else {
    casted <- cleanToNumeric(cleaned1[2:nrow(cleaned1),3:ncol(cleaned1)])
    colnames(casted) <- c("x", "y", as.character(cleaned1[1, 5:ncol(cleaned1)]))
  }
  
  print("removing rows with invalid coordinates")
  casted <- casted[-which(rowSums(is.na(casted[,1:2])) > 0),]
  
  print("removing nan")
  mynatable <- removenan(casted)
  print("imputing na")
  mynatable[,3:ncol(mynatable)] <- halfNA(mynatable[,3:ncol(mynatable)])
  mynatable <- removenan(mynatable)
  
  return(mynatable)
  
  
  
}


#Shows behaviour of a random location, in a specified hour during days, to see difference
#between ferial and festive data
visFest <- function(data, row = sample(1:nrow(data), 1), saturday = 1, hour = 0, period = 24,
                    startCol = 3, endCol = ncol(data)) {
  
  hour <- hour %% period
  svec <- seq(startCol + hour, endCol, period)
  festive <- c(seq(saturday, length(svec), 7), seq(saturday + 1, length(svec), 7))
  
  plot(data[row, seq(startCol + hour, endCol, period)], type = "o", col = "blue", pch = 16,
       main = paste("Daily value, hour ", hour))
  points(festive, data[row, svec[festive]], col = "red", pch = 16)
  
  m <- mean(data[row, svec[festive]])
  s <- sd(data[row,svec[festive]])
  
  abline(h = m, col = "red")
  abline(h = m +s, col = "purple")
  abline(h = m - s, col = "purple")
  
  m2 <- mean(data[row, svec[-festive]])
  s2 <- sd(data[row, svec[-festive]])
  
  abline(h = m2, col = "blue")
  abline(h = m2 + s2, col = "turquoise1")
  abline(h = m2 - s2, col = "turquoise1")
  
}


#Removes rows without valid numbers
#Assumes first two rows are coordinates and does not consider their values
removenan <- function(x) {
  
  #x <- sapply(x, function(y) as.numeric(y))
  
  outRow <- NULL
  for(i in 1:nrow(x)) {
    if(!any(!is.nan(x[i,3:ncol(x)]))) {
      outRow <- c(outRow, i)
    }
    
  }
  #remove fully nan nodes
  if(!is.null(outRow)) {
    
    x <- x[-outRow,]
    
    
  }
  return(x)
}


#Returns best anomalies
topk <- function(results, k, byVotes = FALSE, visual = F, coord = NULL) {
  
  votes <- results$result
  cluster_l <- results$clusters
  
  topkmeasure <- vector("numeric", length = k)
  
  topkelementsT <- vector("numeric", length = k)
  topkelementsI <- vector("numeric", length = k)
  
  if(byVotes) {
    #byvotes
    for(time in 1:length(cluster_l)) {
      for(clust_i in 1:length(cluster_l[[time]])) {
        
        
        clust <- cluster_l[[time]][[clust_i]]
        if(!is.null(clust)){
          voting <- votes[clust$nodes, time]
          mvote <- mean(voting)
          
          if(mvote > topkmeasure[1]) {
            #better priority
            topkmeasure[1] <- mvote
            topkelementsT[1] <- time
            topkelementsI[1] <- clust_i
            
            #reorder
            sorted <- sort(topkmeasure, index.return = T)
            topkelementsT <- topkelementsT[sorted$ix]
            topkelementsI <- topkelementsI[sorted$ix]
            topkmeasure <- sorted$x
          }
        } 
      }
    }
  }
  else {
    #by cluster priority
    for(time in 1:length(cluster_l)) {
      for(clust_i in 1:length(cluster_l[[time]])) {
        
        clust <- cluster_l[[time]][[clust_i]]
        if(!is.null(clust)){
          if(clust$priority > topkmeasure[1]) {
            #better priority
            topkmeasure[1] <- clust$priority
            topkelementsT[1] <- time
            topkelementsI[1] <- clust_i
            
            #reorder
            sorted <- sort(topkmeasure, index.return = T)
            topkelementsT <- topkelementsT[sorted$ix]
            topkelementsI <- topkelementsI[sorted$ix]
            topkmeasure <- sorted$x
          }
        }
      } 
    }
  }
  
  result <- NULL
  for(i in 1:k) {
    
    result <- c(result, list(list(cluster = cluster_l[[topkelementsT[i]]][[topkelementsI[i]]],
                                  measure = topkmeasure[i],
                                  time = topkelementsT[i])))
    
    if(visual) {
      myTimeImage(coord, votes[,topkelementsT[i]], min = 0, max = 40,
                  title = paste(k-i+1,"out of ",k,", time ",topkelementsT[i]))
      myTimeImage.add_cluster(coord, cluster_l[[topkelementsT[i]]][[topkelementsI[i]]],
                              "green")
      #     for(other in 1:length(cluster_l[[topkelementsT[i]]])) {
      #       if(other != topkelementsI[i]) {
      #         myTimeImage.add_cluster(coord, cluster_l[[topkelementsT[i]]][[other]],
      #                                 "green")
      #       }
      #     }
      
    }
  }
  return(result)
}


#Imputes NA with nearest value
#Generates warnings when NA are at the start or end of a series
#Does not impute those NA since they have no meaning
halfNA <- function(x) {
  
  NA_ind <- which(is.na(x))
  for(ind in NA_ind) {
    row <- ind %% nrow(x)
    if(row == 0) {
      row <- nrow(x)
    }
    col <- ceiling(ind / nrow(x))
    
    valid <- which(!is.na(x[row,]))
    
    minBindex <- max(which(valid < col))
    maxAindex <- min(which(valid > col))
    
    if(abs(maxAindex) != Inf && abs(minBindex) != Inf) {
      #don't impute
      
      minBefore <- x[row, valid[max(which(valid < col))]]
      maxAfter <- x[row, valid[min(which(valid > col))]]
      
      beforeDist <- abs(col - valid[max(which(valid < col))])
      afterDist <- abs(col - valid[min(which(valid > col))])
      
      
      x[row,col] <- minBefore + beforeDist * ((maxAfter - minBefore) / (beforeDist + afterDist))
      
    }
    
  }
  return(x)
  
}

#Same as halfNA, with NaN
halfNaN <- function(x) {
  
  NA_ind <- which(is.nan(x))
  for(ind in NA_ind) {
    row <- ind %% nrow(x)
    if(row == 0) {
      row <- nrow(x)
    }
    col <- ceiling(ind / nrow(x))
    
    valid <- which(!is.nan(x[row,]))
    
    minBindex <- max(which(valid < col))
    maxAindex <- min(which(valid > col))
    
    if(abs(maxAindex) != Inf && abs(minBindex) != Inf) {
      #don't impute
      
      minBefore <- x[row, valid[max(which(valid < col))]]
      maxAfter <- x[row, valid[min(which(valid > col))]]
      
      beforeDist <- abs(col - valid[max(which(valid < col))])
      afterDist <- abs(col - valid[min(which(valid > col))])
      
      
      x[row,col] <- minBefore + beforeDist * ((maxAfter - minBefore) / (beforeDist + afterDist))
      
    }
    
  }
  return(x)
  
}



myTimeImage <- function(coord, vals, time, min = NULL, max = NULL,title = NULL, cluster = NULL) {
  
  vals <- vals[,time]
  
  if(is.null(min)){
    min <- min(vals)
  }
  if(is.null(max)) {
    max <- max(vals)
  }
  
  if(is.null(title)) {
    
    title <- paste("Anomalies in time ",time)
  }
  
  vals <- rowMeans(as.matrix(vals), na.rm = T)
  
  mypal <- colorRampPalette( c( "blue", "red" ) )( 5 )
  
  cols <- map2color(vals, mypal, limits = c(min, max))
  
  plot(coord[,1], coord[,2], col = cols, main = title, pch = 19)
  
  if(!is.null(cluster)) {
    if(!is.null(cluster[[time]])) {
      
      
      colList <- c("green", "pink", "purple", "red", "yellow")
      colInd <- 1
      
      for(clu in cluster[[time]]) {
        
        myTimeImage.add_cluster(coord_F, clu, color = colList[colInd])
        colInd <- (colInd) %% length(colList) + 1
      }
        
      
    }
    
    
  }
  
  
  
}

#example
#myTimeImage.add_cluster(coord_F, TEST_6$clusters[[26]][[1]], col = "green")
#
# adds to image first cluster (1) in time 26, with green colors
myTimeImage.add_cluster <- function(coord, cluster, color) {
  
  clust <- cluster$nodes
  points(coord[clust, 1], coord[clust, 2], pch = 19, col = color)
}

slidingWindowsImage <- function(coord, vals, window, min = NULL, max = NULL, overlap = FALSE, save = NULL) {
  
  if(is.null(min)){
    min <- min(vals)
  }
  if(is.null(max)) {
    max <- max(vals)
  }
  
  
  valMat <- as.matrix(vals)
  timeslots <- ncol(valMat)
  
  if(timeslots < window || window < 1) {
    return (NULL)
  }
  
  
  
  if(overlap == TRUE) {
    for(i in (1: (timeslots - window + 1))) {
      if(!is.null(save)) {
        jpeg(file = paste(save,"Image ",i,".jpg"))
      }
      myTimeImage(coord, valMat[,(i:(i + window - 1))], min, max, 
                  title = paste("Mapping values, hours ",i," to ",i+window-1))
      
      
      if(!is.null(save)) {
        dev.off()
      }
    }
  }
  else {
    
    for(i in (0:(floor(timeslots / window) - 1)) ) {
      if(!is.null(save)) {
        jpeg(file = paste(save,"Image ",i+1,".jpg"))
      }
      myTimeImage(coord, valMat[,((i*window + 1):( (i + 1)*window))], min, max, 
                  title = paste("Mapping values, hours ",(i*window + 1)," to ",(i + 1)*window))
      if(!is.null(save)) {
        dev.off()
      }
    }
    
  }
  
  
}

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


# ----- Define a function for plotting a matrix ----- #
#from phaget4.org
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #



densityPlot <- function(data, period = NULL, start = 1) {
  
  if(is.null(period)) {
    vec <- as.numeric(data)
    
  }
  else {
    vec <- as.numeric(data[seq(start, length(data), period)])
    
  }
  
  plot(density(vec, na.rm = T), main = paste("Periodicity ", period, " slot ",start))
  
}









###### Cluster validation ###############





anomaly_cluster.distance5 <- function(anomaly_collection) {
  
  library(dtw)
  
  #find out max and min length of anomaly
  min_length <- Inf
  max_length <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      min_length <- min(min_length, nrow(anomaly))
      max_length <- max(max_length, nrow(anomaly))
      
    }
  }
  
  
  anomalies_by_length <- rep(list(NULL), (max_length - min_length + 1))
  
  feat_num <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      feat_num <- max(feat_num, ncol(anomaly))
      
      alloc <- nrow(anomaly) - min_length + 1
      
      anomalies_by_length[[alloc]] <- c(anomalies_by_length[[alloc]], 
                                        list(anomaly))
      
    }
  }
  
  #if number of anomalies is very high, randomly sample 100 out of them, and calc distances
  
  distance_matrices <- NULL
  max_anoms <- 100
  
  for(i in (1:length(anomalies_by_length))) {
    
    if(!is.null(anomalies_by_length[[i]])) {
      
      if(length(anomalies_by_length[[i]]) > max_anoms) {
        
        anomalies_by_length[[i]] <- anomalies_by_length[[i]][sample(length(anomalies_by_length[[i]]), max_anoms)]
        
      }
      
      feat_distances <- NULL
      all_feat_mat <- NULL
      #build for each feature the distance matrix
      
      for(featN in 1:feat_num) {
        feat_mat <- NULL
        for(anom in anomalies_by_length[[i]]) {
          
          if(ncol(anom) < featN) {
            print("debug")
          }
          
          feat_mat <- rbind(feat_mat, matrix(anom[,featN], nrow = 1))
          
        }
        
        standard_feat <- ( feat_mat - min(feat_mat)) / abs(max(feat_mat) - min(feat_mat))
        
        
        all_feat_mat <- cbind(all_feat_mat, standard_feat)
      }  
      
      
      newDist <- proxy::dist(all_feat_mat, all_feat_mat, method = "correlation")
      
      
      distance_matrices <- c(distance_matrices, list(newDist))
    }
    else {
      distance_matrices <- c(distance_matrices, list(NULL))
    }
  }  
  
  return(list(distances = distance_matrices,
              samples = anomalies_by_length))
  
}


anomaly_cluster.distance4 <- function(anomaly_collection) {
  
  library(dtw)
  
  #find out max and min length of anomaly
  min_length <- Inf
  max_length <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      min_length <- min(min_length, nrow(anomaly))
      max_length <- max(max_length, nrow(anomaly))
      
    }
  }
  
  
  anomalies_by_length <- rep(list(NULL), (max_length - min_length + 1))
  
  feat_num <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      feat_num <- max(feat_num, ncol(anomaly))
      
      alloc <- nrow(anomaly) - min_length + 1
      
      anomalies_by_length[[alloc]] <- c(anomalies_by_length[[alloc]], 
                                        list(anomaly))
      
    }
  }
  
  #if number of anomalies is very high, randomly sample 100 out of them, and calc distances
  
  distance_matrices <- NULL
  max_anoms <- 100
  
  for(i in (1:length(anomalies_by_length))) {
    
    if(!is.null(anomalies_by_length[[i]])) {
      
      if(length(anomalies_by_length[[i]]) > max_anoms) {
        
        anomalies_by_length[[i]] <- anomalies_by_length[[i]][sample(length(anomalies_by_length[[i]]), max_anoms)]
        
      }
      
      feat_distances <- NULL
      #build for each feature the distance matrix
      
      for(featN in 1:feat_num) {
        feat_mat <- NULL
        for(anom in anomalies_by_length[[i]]) {

          
          feat_mat <- rbind(feat_mat, matrix(anom[,featN], nrow = 1))
          
        }
        
        #before distance, standardize feature
        standard_feat <- ( feat_mat - min(feat_mat)) / abs(max(feat_mat) - min(feat_mat))
        
        newDist <- proxy::dist(standard_feat, standard_feat,method = "correlation")
        
        feat_distances <- cbind(feat_distances, newDist)
        
      }
      
      distance_matrices <- c(distance_matrices, list(feat_distances))
    }
    else {
      distance_matrices <- c(distance_matrices, list(NULL))
    }
  }  
  
  return(list(distances = distance_matrices,
              samples = anomalies_by_length))
  
}



anomaly_cluster.distance3 <- function(anomaly_collection) {
  
  library(dtw)
  
  #find out max and min length of anomaly
  min_length <- Inf
  max_length <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      min_length <- min(min_length, nrow(anomaly))
      max_length <- max(max_length, nrow(anomaly))
      
    }
  }
  
  
  anomalies_by_length <- rep(list(NULL), (max_length - min_length + 1))
  
  feat_num <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      feat_num <- max(feat_num, ncol(anomaly))
      
      alloc <- nrow(anomaly) - min_length + 1
      
      anomalies_by_length[[alloc]] <- c(anomalies_by_length[[alloc]], 
                                        list(anomaly))
      
    }
  }
  
  #if number of anomalies is very high, randomly sample 100 out of them, and calc distances
  
  distance_matrices <- NULL
  max_anoms <- 100
  
  for(i in (1:length(anomalies_by_length))) {
    
    if(!is.null(anomalies_by_length[[i]])) {
      
      if(length(anomalies_by_length[[i]]) > max_anoms) {
        
        anomalies_by_length[[i]] <- anomalies_by_length[[i]][sample(length(anomalies_by_length[[i]]), max_anoms)]
        
      }
      
      feat_distances <- NULL
      all_feat_mat <- NULL
      #build for each feature the distance matrix
      
      for(featN in 1:feat_num) {
        feat_mat <- NULL
        for(anom in anomalies_by_length[[i]]) {
          
          if(ncol(anom) < featN) {
            print("debug")
          }
          
          feat_mat <- rbind(feat_mat, matrix(anom[,featN], nrow = 1))
          
        }
        
        standard_feat <- ( feat_mat - min(feat_mat)) / abs(max(feat_mat) - min(feat_mat))
        
        
        all_feat_mat <- cbind(all_feat_mat, standard_feat)
      }  
      

        newDist <- proxy::dist(all_feat_mat, all_feat_mat, method = "DTW")
        
      
      distance_matrices <- c(distance_matrices, list(newDist))
    }
    else {
      distance_matrices <- c(distance_matrices, list(NULL))
    }
  }  
  
  return(list(distances = distance_matrices,
              samples = anomalies_by_length))
  
}


anomaly_cluster.distance2 <- function(anomaly_collection) {
  
  library(dtw)
  
  #find out max and min length of anomaly
  min_length <- Inf
  max_length <- 0
  
  for(node_list in anomaly_collection) {
      for(anomaly in node_list) {
        
        min_length <- min(min_length, nrow(anomaly))
        max_length <- max(max_length, nrow(anomaly))
        
      }
  }
  
  
  anomalies_by_length <- rep(list(NULL), (max_length - min_length + 1))
  
  feat_num <- 0
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      feat_num <- max(feat_num, ncol(anomaly))
      
      alloc <- nrow(anomaly) - min_length + 1
      
      anomalies_by_length[[alloc]] <- c(anomalies_by_length[[alloc]], 
                                        list(anomaly))
      
    }
  }
  
  #if number of anomalies is very high, randomly sample 100 out of them, and calc distances
  
  distance_matrices <- NULL
  max_anoms <- 100
  
  for(i in (1:length(anomalies_by_length))) {
    
    if(!is.null(anomalies_by_length[[i]])) {
      
      if(length(anomalies_by_length[[i]]) > max_anoms) {
        
        anomalies_by_length[[i]] <- anomalies_by_length[[i]][sample(length(anomalies_by_length[[i]]), max_anoms)]
        
      }
      
      feat_distances <- NULL
      #build for each feature the distance matrix
      
      for(featN in 1:feat_num) {
        feat_mat <- NULL
        for(anom in anomalies_by_length[[i]]) {
          
          if(ncol(anom) < featN) {
            print("debug")
          }
          
          feat_mat <- rbind(feat_mat, matrix(anom[,featN], nrow = 1))
          
        }
        
        #before distance, standardize feature
        standard_feat <- ( feat_mat - min(feat_mat)) / abs(max(feat_mat) - min(feat_mat))
        
        newDist <- proxy::dist(standard_feat, standard_feat,method = "DTW")
            
        feat_distances <- cbind(feat_distances, newDist)
        
      }
      
      distance_matrices <- c(distance_matrices, list(feat_distances))
    }
    else {
      distance_matrices <- c(distance_matrices, list(NULL))
    }
  }  
  
  return(list(distances = distance_matrices,
              samples = anomalies_by_length))
  
}



anomaly_cluster.distance <- function(anomaly_collection) {
  
  library(dtw)
  
  #find out max and min length of anomaly
  min_length <- Inf
  max_length <- 0
  
    for(node_list in anomaly_collection) {
      for(anomaly in node_list) {
      
        min_length <- min(min_length, length(anomaly))
        max_length <- max(max_length, length(anomaly))
      
      }
    }
  
  
  anomalies_by_length <- rep(list(NULL), (max_length - min_length + 1))
  
  
  for(node_list in anomaly_collection) {
    for(anomaly in node_list) {
      
      alloc <- length(anomaly) - min_length + 1
      
      anomalies_by_length[[alloc]] <- rbind(anomalies_by_length[[alloc]],
                                            anomaly)
      
    }
  }
  
  #if number of anomalies is very high, randomly sample 100 out of them, and calc distances
  
  distance_matrices <- NULL
  max_rows <- 100
  
  for(i in (1:length(anomalies_by_length))) {
    
    if(!is.null(anomalies_by_length[[i]])) {
      
      if(nrow(anomalies_by_length[[i]]) > max_rows) {
        
        anomalies_by_length[[i]] <- anomalies_by_length[[i]][sample(nrow(anomalies_by_length[[i]]), max_rows),]
        
      }
      
      distance_matrices <- c(distance_matrices, list(proxy::dist(anomalies_by_length[[i]],
                                                                 anomalies_by_length[[i]],
                                                                 method = "DTW")))
      
    }
    else {
      distance_matrices <- c(distance_matrices, list(NULL))
    }
  }  
  
  return(list(distances = distance_matrices,
              samples = anomalies_by_length))
  
}

anomaly_cluster.clustering <- function(distance_results, k) {
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    if(!is.null(distances[[i]]) && nrow(samples[[i]]) > k) {
      
      
      clust <- c(clust, list(kmeans(distances[[i]], centers = k)))
      
    } 
    else {
      clust <- c(clust, list(samples[[i]]))
    }
  }
  
  return(clust)
}


anomaly_cluster.bestClustering3 <- function(distance_results) {
  library(fpc)
  
  test <- NULL
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    
    maxK <- min(30, length(samples[[i]]) - 1)
    
    if(!is.null(distances[[i]]) && length(distances[[i]]) > 1) {
      best_clust <- pamk(distances[[i]], krange = (2:maxK))
      
  
      clust <- c(clust, list(best_clust))
      
    }
    else {
      clust <- c(clust, list(samples[[i]]))
    }
    
  }
  
  plot(2:length(test), test[-1], type = "o")
  
  return(clust)
}

anomaly_cluster.bestClustering2 <- function(distance_results) {
  library(cluster)
  
  test <- NULL
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    
    maxK <- min(30, length(samples[[i]]) - 1)
    
    if(!is.null(distances[[i]]) && length(distances[[i]]) > 1) {
      best_clust <- NULL
      for(k in 1:maxK) {
        tent_clust <- pam(distances[[i]], k)
        test <- c(test, abs(mean(tent_clust$tot.withinss / tent_clust$betweenss)))
        
        if(is.null(best_clust)){
          best_clust <- tent_clust
        }
        else if(abs(mean(tent_clust$tot.withinss / tent_clust$betweenss)) < abs(mean(best_clust$tot.withinss / best_clust$betweenss))) {
          
          best_clust <- tent_clust
        }
        
      }
      
      clust <- c(clust, list(best_clust))
      
    }
    else {
      clust <- c(clust, list(samples[[i]]))
    }
    
  }
  
  plot(2:length(test), test[-1], type = "o")
  
  return(clust)
}


anomaly_cluster.clustering2 <- function(distance_results, k) {
  
  library(cluster)
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    if(!is.null(distances[[i]]) && length(samples[[i]]) > k) {
      
      
      clust <- c(clust, list(pam(distances[[i]], k)))
      
    } 
    else {
      clust <- c(clust, list(samples[[i]]))
    }
  }
  
  return(clust)
}

anomaly_cluster.clustering3 <- function(distance_results) {
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    if(!is.null(distances[[i]])) {
      
      
      clust <- c(clust, list(hclust(distances[[i]])))
    } 
    else {
      clust <- c(clust, list(samples[[i]]))
    }
  }
  
  return(clust)
}

anomaly_cluster.showClust <- function(cluster, visual = TRUE) {
   
    mean <- 0
    
    for(samplesI in 1:length(cluster)) {
        
      mean <- ( (mean*(samplesI - 1) + cluster[[samplesI]]) )/ samplesI
      
    }
    
    if(visual) {
      for(f in 1:ncol(mean)) {
        plot(mean[,f] ,main = paste("Feature ", f), type = "l")
      }
    }
    return(mean)
  
}

anomaly_cluster.show2 <- function(anomalies, feat_numb, k = 3) {
  
  for(feat in feat_numb) {
    
    anom_F.dist <- anomaly_cluster.distance(anomalies[[feat]])
    anom_F.clust <- anomaly_cluster.clustering(anom_F.dist, k)
    
    for(i in (1:k)) {
      
      clustwindows <- anom_F.dist$samples[[1]][which(anom_F.clust[[1]]$cluster == i),]
      plot(colMeans(clustwindows), type = "l", main = paste("Type ",i,", feature ",feat))      
    }
    
    
  }
  
  
}

anomaly_cluster.bestClustering <- function(distance_results) {
  
  test <- NULL
  
  distances <- distance_results$distances
  samples <- distance_results$samples
  
  clust <- NULL
  
  for(i in (1:length(distances))) {
    
    maxK <- min(30, nrow(samples[[i]]) - 1)
    
    if(!is.null(distances[[i]]) && length(distances[[i]]) > 1) {
      best_clust <- NULL
      for(k in 1:maxK) {
        tent_clust <- kmeans(distances[[i]], centers = k)
        test <- c(test, abs(mean(tent_clust$tot.withinss / tent_clust$betweenss)))
        
        if(is.null(best_clust)){
          best_clust <- tent_clust
        }
        else if(abs(mean(tent_clust$tot.withinss / tent_clust$betweenss)) < abs(mean(best_clust$tot.withinss / best_clust$betweenss))) {
          
          best_clust <- tent_clust
        }
        
      }
      
      clust <- c(clust, list(best_clust))
      
    }
    else {
      clust <- c(clust, list(samples[[i]]))
    }
    
  }
  
  plot(2:length(test), test[-1], type = "o")
  
  return(clust)
}

anomaly_cluster.show <- function(clust, samples, clk = NULL, numb = NULL) {
  
  if(is.null(numb)) {
    #number of width to be considered
    numb <- 1:length(clust)
  }
  if(is.null(clk)) {
    #number of cluster types to consider
    clk <- 1:length(unique(clust[[1]]$cluster))
    
  }
  
  
  plotRows <- length(clk)
  plotCols <- length(numb)
  
  
  old.par <- par(mfrow = c(plotRows,plotCols), mai = c(0.6,1,0.5,0.5))
  
  for(i in numb) {
    for(j in clk) {
      
      indices <- which(clust[[i]]$cluster == j)
      samp <- samples[[i]][indices,]
      main <- paste("Anomaly type ",j, ", length ",i)
      
      if(length(indices) == 1) {
        plot(samp, type = "l", main = main)
      }
      else {
        plot(colMeans(samp), type = "l", main = main)
      }
    }
  }
  
  par(old.par)
}




###################### K-MEANS #########


#myData: a matrix, each row corresponding to a window of a different time series
#function returns the number of clusters that best fit the data
#algorithm used is k-medoid, a robust variant of k-means
#Works with windows of width 1, for a single time slot analysis
bestClusterNumber <- function(myData, verbose = FALSE) {
  
  library(dtw) #for dist
  
  library(fpc) # for pamk
  
  if(verbose == TRUE) {
    print("Computing distances...")
  }
  distMatrix <- dist(myData, method = "DTW")
  if(verbose == TRUE) {
    print("Done computing distances")
  }
  
  myFit <- pamk(distMatrix, krange = (1:(nrow(distMatrix) - 1)))
  
  if(verbose == TRUE) {
    print(paste("Pamk says... ",myFit$nc,"!"))
    if(myFit$nc > 2) {
      print(paste("Pamk thus spach"))
      print(myFit[[1]]$clustering)
    }
    plot(myFit$pamobject, main = "Pamk!!!!!")
    
    
    
  }
  
  return(myFit$nc)
  
  
  return(nclus)
}

FOM <- function(fit) {
  
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  
  AIC <- (D + log(n)*m*k)
  
  return(AIC)
  
}

testBest <- function(myData, window) {
  
  verb = TRUE
  
  
  result <- NULL
  
  if(window > 1) {
    for(i in (1:window - 1)) {
      
      result <- c(result, 0)
      
      
    }
    
  }
  
  for(i in c(window:ncol(myData))) {
    
    print(paste("Time ",i - window + 1," to time ",i))
    result <- c(result,
                bestClusterNumber(myData[,c((i - window + 1):i)], 
                                  verbose = verb))
    
  }
  
  return(result)
}

####### FLUFF #######

pwr <- function(x, exp) {
  
  return (sign(x) * abs(x) ^ exp)
  
  
}

showFits <- function(meth_fit, node, kpi, method = NULL) {
  
  if(is.null(method)){
    method <- length(meth_fit[[node]][[kpi]])
    
  }
  
  for(i in method) {
    
    dense <- density(meth_fit[[node]][[kpi]][[i]])
    plot(dense, main = paste("Density, node ",node,", kpi ",kpi,", method ",i))
    
  }
  
  
}

SMOTE <- function() {
  
}

anomaly_windows <- function(data, anomaly_times, window) {
  
  
  #find anomaly sets
  current_motif <- NULL
  motif_set <- NULL
  
  for(cand in (anomaly_times)) {
      
    for(feat in data) {
      
      current_motif <- c(data[[feat]])
        
      
    }
    
  
  }
    
  
  
}

anomaly_window_search <- function(data, anomaly_times, min_window) {
  
  #find anomaly sets
  
  current_motif <- NULL
  motif_set <- NULL
  
  for(cand in (anomaly_times)) {
    
    #find motif representing found anomaly
    current_motif <- c(current_motif, data[cand])
    
    if( !((cand + 1) %in% anomaly_times)) {
      #stored all anomalous moments
      
      before <- 0
      after <- 0
      noMore <- F
      init_length <- length(current_motif)
      
      while(length(current_motif) < min_window && noMore == FALSE) {
        #pad data
        
        if(before <= after && (cand - init_length - before) >= 1 ) {
          current_motif <- c(data[cand - init_length - before], current_motif)
          before <- before + 1
        }
        else if(cand + after + 1 <= length(data))   {
          current_motif <- c(current_motif, data[cand + after + 1])
          after <- after + 1
        }
        else if((cand - init_length - before) >= 1 ) {
          current_motif <- c(data[cand - init_length - before], current_motif)
          before <- before + 1
        }
        else {
          noMore <- TRUE
        }
      }
      
      motif_set <- c(motif_set, list(current_motif))
      current_motif <- NULL
      
    }
  }
  
  return(motif_set)
  
}