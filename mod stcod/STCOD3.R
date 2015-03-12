
######## Main Algorithm ##########################################

#cellData: list of n kpi, each containing a matrix for all element
#neigh_number: includes self in count of neighborhood
STCOD_Cluster2 <- function(cellData, trainingData = NULL, coordinates, period = 24, watch_node = 5, neigh_number = 20, visual = FALSE, 
                           methodAct = NULL, motif_discovery = FALSE, NA_strat = "local_mean",
                           anomaly_window = 4, anomaly_window_collection = TRUE, festive) {
  
  state_mem <- 5
  
  confidence_threshold <- 0.997
  
  memory_periods <- 3
  history_periods <- 3
  
  min_difference <- 0.0001
  
  kpi_number <- length(cellData)
  node_number <- nrow(cellData[[1]])
  final_timeslot <- ncol(cellData[[1]])
  
  
  
  method_number = 5
  selectMethod <- vector("logical", length = method_number)
  
  #ACTIVATE\DEACTIVATE METHOD
  if(!is.null(methodAct)) {
    selectMethod <- methodAct
  }
  else {
    selectMethod <- rep(TRUE, method_number)
  }
  
  
  startSeasonCheck <- period + 1
  startTrendCheck <- startSeasonCheck
  startHistoryCheck <- (history_periods * period) + 1
  startMemoryCheck <- (memory_periods * period) + 1
  
  #indicates when to stop training period for anomaly forecasting (motif discovery will be done the timeslot preceding this)
  anomaly_testing_time <- (9 * period) + 1
  
  distance_threshold <- 1
  descent_threshold <- 1
  threshold_confidence <- 0.63
  
  #indicates quality of method
  method_weight <- rep(list(matrix(1, nrow = method_number, ncol = kpi_number)), node_number)
  
  anomaly_significance <- 0.997
  max_priority <- 0.999999
  
  
  
  #following data is stored in list instead of matrix to simulate cell distributed behavior
  #element in a list[[1]] represents knowledge of node 1
  #this is redundant in this algorithm, as cluster_season_sample and prev_state should be the same
  
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
  cell_memories <- rep(list(matrix(rep(0, kpi_number *(period*memory_periods)), 
                                   nrow = kpi_number, 
                                   ncol = (period * memory_periods))), node_number)
  
  #historical values
  cell_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  
  cell_descent_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_descent_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_distance_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_distance_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  #state memory of cell and neighbors
  #prev_state[[node]][1] indicates current state (after phase 1), [2] indicates state in previous time etc...
  prev_state <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
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
  node_motifs <- rep(list(NULL), node_number)
  anomaly_motifs <- rep(list(NULL), node_number)
  anomaly_collection <- rep(list(rep(list(NULL), node_number)), kpi_number)
  anomaly_times <- rep(list(NULL), node_number)
  priority_history <- matrix(,nrow = node_number, ncol = (anomaly_testing_time - 1))
  
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
    
    cell_descent_historical_mean <- trainingData$descent_mean
    cell_descent_historical_sd <- trainingData$descent_sd
    cell_distance_historical_mean <- trainingData$distance_mean
    cell_distance_historical_sd <- trainingData$distance_sd
    
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
  
  #find bestK to partition a time series
  
  if(skipTraining == T) {
    
    startSeasonCheck <- 1
    startHistoryCheck <- 1
    startTrendCheck <- 2
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
      
      
      if(selectMethod[1] == TRUE) {
        season_fit <- vector("numeric", length = kpi_number)
      }
      if(selectMethod[2] == TRUE) {
        distance_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[3] == TRUE) {
        trend_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[4] == TRUE) {
        previous_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[5] == TRUE) {
        history_fit <- vector("numeric", length = kpi_number)
        
      }
      
      #for each feature in the node, calculate various measures of unlikelihood
      for(kpi in (1:kpi_number)) {
        
        curr_values[kpi] <- cellData[[kpi]][node, time]
        
        if(festivity) {
          curr_seasonality <- cell_F_seasonality[[node]][kpi,]
          
        }
        else {
          curr_seasonality <- cell_seasonality[[node]][kpi,]
          
        }
        
        if(time >= startSeasonCheck) {
          
          if(selectMethod[1] == TRUE) {
            season_fit[kpi] <- featureSeasonFitness(curr_seasonality,
                                                    curr_values[kpi],
                                                    curr_period, replicas = 3, visual = visual)
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
          
          
          if(skipTraining == FALSE) {
            
            
            last_distance_mean <- cell_distance_historical_mean[[node]][kpi]
            
            cell_distance_historical_mean[[node]][kpi] <- (cell_distance_historical_mean[[node]][kpi] * (time - 1) +
                                                             dist) / time
            cell_distance_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_distance_historical_sd[[node]][kpi] ^2) ) 
                                                              + (dist - last_distance_mean) * (dist - cell_distance_historical_mean[[node]][kpi])) / time)
            
            if(is.nan(cell_distance_historical_sd[[node]][kpi])) {
              cell_distance_historical_sd[[node]][kpi] <- 0
            }
          }
          
          
          
          
        }
        
        if(time >= startTrendCheck) {
          
          
          prev_val_std <- (cell_memories[[node]][kpi,ncol(cell_memories[[node]])] - season_avg) / season_rng 
          prev_season_std <- (curr_seasonality[last_period] - season_avg) / season_rng 
          
          residual <- (curr_val_std - prev_val_std) - (curr_season_std - prev_season_std)
          
          if(skipTraining == FALSE) {
            
            last_descent_mean <- cell_descent_historical_mean[[node]][kpi]
            
            cell_descent_historical_mean[[node]][kpi] <- (cell_descent_historical_mean[[node]][kpi] * (time - 1) +
                                                            residual) / time
            cell_descent_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_descent_historical_sd[[node]][kpi] ^2) ) 
                                                             + (residual - last_descent_mean) * (residual - cell_descent_historical_mean[[node]][kpi])) / time)
            if(is.nan(cell_descent_historical_sd[[node]][kpi])) {
              cell_descent_historical_sd[[node]][kpi] <- 0
            }
            
          }
          
          if(selectMethod[3] == TRUE) {
            
            trend_fit[kpi] <- featureTrendThresholdPVal(curr_season_std,
                                                        prev_season_std,
                                                        curr_val_std,
                                                        prev_val_std,
                                                        cell_descent_historical_sd[[node]][kpi],
                                                        threshold_confidence, visual)
            
          }
          
        }
        if(time >= startMemoryCheck) {
          
          if(selectMethod[4] == TRUE) {
            previous_fit[kpi] <- featureLastValuesFitness(cell_memories[[node]][kpi,],
                                                          curr_values[kpi], visual)          
          }
        }
        
        if(time >= startHistoryCheck) {
          
          if(cell_historical_sd[[node]][kpi] < min_difference) {
            cell_historical_sd[[node]][kpi] <- min_difference
          }
          
          if(selectMethod[2] == TRUE) {
            distance_fit[kpi] <- featureDistanceThresholdPVal(curr_season_std,
                                                              curr_val_std,
                                                              cell_distance_historical_sd[[node]][kpi],
                                                              threshold_confidence, visual)
          }
          
          
          
          if(selectMethod[5] == TRUE) {
            if(festivity) {
              history_fit[kpi] <- featureHistoryFitness(cell_F_historical_mean[[node]][kpi],
                                                        cell_F_historical_sd[[node]][kpi],
                                                        curr_values[kpi], visual)
              
            }
            else {
              history_fit[kpi] <- featureHistoryFitness(cell_historical_mean[[node]][kpi],
                                                        cell_historical_sd[[node]][kpi],
                                                        curr_values[kpi], visual)
              
            }
          }
          
        }
        else {
          
          
          #update historical mean and deviation for kpi
          last_mean <- cell_historical_mean[[node]][kpi]
          
          cell_historical_mean[[node]][kpi] <- (cell_historical_mean[[node]][kpi] * (time - 1) +
                                                  curr_values[kpi]) / time
          cell_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_historical_sd[[node]][kpi] ^2) ) 
                                                   + (curr_values[kpi] - last_mean) * (curr_values[kpi] - cell_historical_mean[[node]][kpi])) / time)
          
          
        }
        
        
      }
      
      
      
      #if no comparison value exists, skip any check, consider node OK (learning phase)
      if(time < min(startSeasonCheck, startHistoryCheck, startMemoryCheck)) {
        
        #update current state, shift previous for this node in its list
        prev_state[[node]][1,] <- shiftRight(prev_state[[node]][1,] , 1)
        prev_state[[node]][1,1] <- "OK"
        
        if(skipTraining == FALSE) {
          
          cell_seasonality[[node]] <- updateSeasonality(cell_seasonality[[node]], curr_values, curr_period,
                                                        curr_status = prev_state[[node]][1,1], noSeasonPresent = TRUE)
          
        }
        
        
        node_priority[[node]][1] <- 0
        
        prev_priority[[node]][1,] <-  shiftRight(prev_priority[[node]][1,] , 1)
        prev_priority[[node]][1,1] <- node_priority[[node]][1]
        
      }
      #else check values
      else {
        
        methodisvalid <- (!(time < c(startSeasonCheck, startHistoryCheck, startHistoryCheck, startMemoryCheck, startHistoryCheck))) & (selectMethod)
        
        #aggregate by kpi
        univariateEnsemblePriority <- vector("numeric", length = kpi_number)
        
        for(kpi in (1:kpi_number)) {
          
          kpiPVals <- NULL
          
          if(time >= startSeasonCheck) {
            if(selectMethod[1] == TRUE) {
              kpiPVals <- c(kpiPVals,
                            season_fit[kpi])
            }
            if(selectMethod[2] == TRUE && time >= startHistoryCheck) {
              
              kpiPVals <- c(kpiPVals,
                            distance_fit[kpi])
            }            
            if(selectMethod[3] == TRUE && time >= startTrendCheck) {
              
              kpiPVals <- c(kpiPVals,
                            trend_fit[kpi])
            }
          }
          if(time >= startMemoryCheck) {
            if(selectMethod[4] == TRUE) {
              kpiPVals <- c(kpiPVals,
                            previous_fit[kpi])
            }
          }
          if(time >= startHistoryCheck) {
            if(selectMethod[5] == TRUE) {
              kpiPVals <- c(kpiPVals,
                            history_fit[kpi])
            }
          }
          
          
          univariateEnsemblePriority[kpi] <- fisherCombinedPriority(kpiPVals)
          
        }
        
        #aggregate by method
        multivariateMethodPriority <- NULL 
        
        if(time >= startSeasonCheck) {
          if(selectMethod[1] == TRUE) {
            multivariateMethodPriority <- c(fisherCombinedPriority(season_fit))
          }
          if(selectMethod[2] == TRUE && time >= startHistoryCheck) {
            
            multivariateMethodPriority <- c(multivariateMethodPriority,
                                            fisherCombinedPriority(distance_fit))
          }
          if(time >=startTrendCheck) {
            if(selectMethod[3] == TRUE) {
              multivariateMethodPriority <- c(multivariateMethodPriority,
                                              fisherCombinedPriority(trend_fit))
              
            }
          }
        }
        if(time >= startMemoryCheck) {
          if(selectMethod[4] == TRUE) {
            multivariateMethodPriority <- c(multivariateMethodPriority, 
                                            fisherCombinedPriority(previous_fit))
          }
        }
        if(time >= startHistoryCheck) {
          if(selectMethod[5] == TRUE) {
            multivariateMethodPriority <- c(multivariateMethodPriority, 
                                            fisherCombinedPriority(history_fit))
          }
        }
        
        if(time >= 73 && is.nan(multivariateMethodPriority)) {
          print("debug")
        }
        
        if(any(multivariateMethodPriority > confidence_threshold)) {
          
          for( method in (1: length(multivariateMethodPriority))) {
            
            if(multivariateMethodPriority[method] > confidence_threshold) {
              
              #print(paste("Anomaly node ",node," time ", time," signaled by method ", method ," over all kpi"))
              prev_state[[node]][1,] <- shiftRight(prev_state[[node]][1,] , 1)
              prev_state[[node]][1,1] <- "Anomaly"
            }
            
          }
          
          
        }
        else if(any(univariateEnsemblePriority > confidence_threshold)) {
          print(paste("Anomaly node ",node," time ", time," signaled over a kpi by ensemble of methods"))
          prev_state[[node]][1,] <- shiftRight(prev_state[[node]][1,] , 1)
          prev_state[[node]][1,1] <- "Anomaly"
        }
        else {
          prev_state[[node]][1,] <- shiftRight(prev_state[[node]][1,] , 1)
          prev_state[[node]][1,1] <- "OK"
        }
        
        alfa <- confidence_threshold
        
        node_priority[[node]][1] <- priority_function(1 - univariateEnsemblePriority, alfa)
        
        prev_priority[[node]][1,] <-  shiftRight(prev_priority[[node]][1,] , 1)
        prev_priority[[node]][1,1] <- node_priority[[node]][1]
        
        if(node_priority[[node]][1] > anomaly_significance) {
          
          print(paste("Anomalous node ",node,", overall"))
        }
        
        #update seasonality
        if(skipTraining == FALSE) {
          cell_seasonality[[node]] <- updateSeasonality(cell_seasonality[[node]], curr_values, curr_period,
                                                        curr_status = prev_state[[node]][1,1])
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
          
          prev_state[[node]][neigh_ind,] <- shiftRight(prev_state[[node]][neigh_ind,], 1)
          prev_state[[node]][neigh_ind,1] <- prev_state[[true_neighbour_index]][1,1]
          
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
      
      if(time >= max(min(startSeasonCheck, startMemoryCheck, startHistoryCheck), 1 + startTrendCheck + ncol(prev_priority[[node]]))) {
        
        
        cluster_op <- fastGeneralWeightedFisherScan(node_priority[[node]], neigh_weight[[node]], anomaly_significance)
        
        if (cluster_op$vote > anomaly_significance) {
          
          #print(paste("Node: ",node,", anomaly cluster at time ",time))
          
          
          inv_neighs <- cluster_op$nodes
          inv_nodes <- NULL
          
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
      
      
      #SIMILARITY CHECK - WEIGHT UPDATING
      
      if(time >= startTrendCheck + ncol(prev_priority[[node]])) {
        similarity_update_rate <- 0.2
        
        for(neigh_ind in (2:neigh_number)) {
          
          
          similarity <- similarity_check(prev_priority[[node]][1,], prev_priority[[node]][neigh_ind,])
          
          neigh_weight[[node]][neigh_ind] <- ((1 - similarity_update_rate)*neigh_weight[[node]][neigh_ind] 
                                              + similarity_update_rate * similarity)
          
        }
        
        
      }
      
      #ANOMALY FORECASTING
      
      
      if(time < anomaly_testing_time) {
        
        priority_history[node, time] <- node_priority[[node]][1]
        
      }
      
      if(motif_discovery == T && time == anomaly_testing_time - 1 && !is.null(anomaly_times[[node]])) {
        
        priority_history[is.na(priority_history)] <- 0
        
        #if not known, find bestK
        if(is.null(bestK)) {
          
          print("labeling for best k")
          #labeled_set <- psf_label_set(cellData[[1]][neigh_index[[node]], 1:time], 4, 24)
          labeled_set <- psf_label_set(priority_history[neigh_index[[node]], 1:time], 4, 24)
          
          allk <- 0
          maxK <- 0
          for(i in 1:length(labeled_set$levels)) {
            allk <- allk + length(labeled_set$levels[[i]])
            maxK <- max(maxK, length(labeled_set$levels[[i]]))
            
          }
          
          #bestK <- ceiling(allk/length(labeled_set$levels) )
          bestK <- maxK
          
        }
        
        #for now 1 kpi
        #new_motifs <-  motif_anomaly_discovery(cellData[[1]][neigh_index[[node]], 1:time], 
        #                                       anomaly_times[[node]],
        #                                       5, min_window = 10, max_window = 10,
        #                                       min_k = bestK, max_k = bestK)
        
        
        # new_chains <- anomaly_motifs(cellData[[1]][neigh_index[[node]], 1:time], 
        #                               anomaly_times[[node]], min_k = bestK, max_k = bestK)
        new_chains <- anomaly_motifs(priority_history[neigh_index[[node]], 1:time], 
                                     anomaly_times[[node]], min_k = bestK, max_k = bestK)
        
        #         
        #         if(!is.null(new_motifs)) {
        #           
        #           node_motifs[[node]] <- new_motifs
        #         }
        if(!is.null(new_chains$set)) {
          
          anomaly_motifs[[node]] <- new_chains$set        
          bestLevels[[node]] <- new_chains$lev
          
        }
        
      }
      else if(motif_discovery == T && time >= anomaly_testing_time && !is.null(anomaly_motifs[[node]])) {
        #Search through the discovered motifs
        #        forecasts <- motif_search(as.matrix(cellData[[1]][neigh_index[[node]], time]), bestLevels[[node]],
        #                                  anomaly_motifs[[node]], bestK)
        forecasts <- motif_search(as.matrix(node_priority[[node]]), bestLevels[[node]],
                                  anomaly_motifs[[node]], bestK)
        
        if(any(forecasts > 0)) {
          
          times <- which(forecasts > 0)
          probs <- forecasts[times]
          string <- paste(times - 1, "hours with prob", probs,";", collapse = " ")
          print(paste("Node ",node," forecasting anomaly in ", string))
          
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
    
    #     if(any(node_votes > 3)) {
    #       #       
    #       #       library(fpc)
    #       #       pamcl <- pamk(node_votes[node_votes > 3], krange = (1:(length(node_votes[node_votes > 3]) - 1)))
    #       #       
    #       print(which(node_votes > 3))
    #       
    #     }
    
    
    final_result <- cbind(final_result, node_votes)
    
  }
  
  print("end")
  
  if(anomaly_window_collection == TRUE) {
    for(node in 1:node_number) {
      
      anomaly_times <- which(final_result[node,] > 0)
      
      for(feature in 1:length(cellData)) {
        found_anomalies <- anomaly_window_search(as.numeric(cellData[[feature]][node,]), anomaly_times, anomaly_window)
        if(!is.null(found_anomalies)) {
          anomaly_collection[[feature]][[node]] <- found_anomalies
        }
      }
    }
  }
  
  return(list(result = final_result,
              anomalies = anomaly_collection,
              clusters = anomaly_clusters))
  
}




######## Training Algorithm ##############



STCOD_Cluster2.training <- function(cellData, festive, coordinates, period = 24, neigh_number = 10,
                                    methodAct = NULL, NA_strat = "local_mean") {
  #for kernel cdf
  library(sROC)
  
  state_mem <- 5
  
  confidence_threshold <- 0.997
  
  
  memory_periods <- 3
  history_periods <- 3
  
  min_difference <- 0.0001
  
  kpi_number <- length(cellData)
  node_number <- nrow(cellData[[1]])
  final_timeslot <- ncol(cellData[[1]])
  
  
  method_number = 5
  selectMethod <- vector("logical", length = method_number)
  
  #ACTIVATE\DEACTIVATE METHOD
  if(!is.null(methodAct)) {
    selectMethod <- methodAct
  }
  else {
    selectMethod <- rep(FALSE, method_number)
  }
  
  
  startSeasonCheck <- period + 1
  startHistoryCheck <- (history_periods * period) + 1
  startMemoryCheck <- (memory_periods * period) + 1
  
  #indicates when to stop training period for anomaly forecasting (motif discovery will be done the timeslot preceding this)
  anomaly_testing_time <- (9 * period) + 1
  
  distance_threshold <- 1
  descent_threshold <- 1
  threshold_confidence <- 0.842
  
  #indicates quality of method
  method_weight <- rep(list(matrix(1, nrow = method_number, ncol = kpi_number)), node_number)
  
  anomaly_significance <- 0.9545
  max_priority <- 0.999999
  
  
  
  #following data is stored in list instead of matrix to simulate cell distributed behavior
  #element in a list[[1]] represents knowledge of node 1
  #this is redundant in this algorithm, as cluster_season_sample and prev_state should be the same
  
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
  cell_memories <- rep(list(matrix(rep(0, kpi_number *(period*memory_periods)), 
                                   nrow = kpi_number, 
                                   ncol = (period * memory_periods))), node_number)
  
  #historical values
  cell_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_F_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_F_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  
  cell_descent_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_descent_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_distance_historical_mean <- rep(list(vector("numeric", length = kpi_number)), node_number)
  cell_distance_historical_sd <- rep(list(vector("numeric", length = kpi_number)), node_number)
  
  cell_cdf <- rep(list(rep(list(NULL), kpi_number)), node_number)
  
  
  #state memory of cell and neighbors
  #prev_state[[node]][1] indicates current state (after phase 1), [2] indicates state in previous time etc...
  prev_state <- rep(list(matrix(, nrow = neigh_number, ncol = state_mem)), node_number)
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
  node_motifs <- rep(list(NULL), node_number)
  anomaly_motifs <- rep(list(NULL), node_number)
  anomaly_collection <- rep(list(NULL), node_number)
  anomaly_times <- rep(list(NULL), node_number)
  priority_history <- matrix(,nrow = node_number, ncol = (anomaly_testing_time - 1))
  
  #because this implementation is not distributed, to avoid computing n*k times a bestK clustering
  #only the first node does, and then it's assumed valid for all
  
  bestK <- NULL
  bestLevels <- rep(list(NULL), node_number)
  
  
  ###################Algorithm start
  
  
  
  
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
      #       outRow <- NULL
      #       for(i in 1:nrow(cellData[[feature]])) {
      #         if(!any(!is.nan(cellData[[feature]][i,3:ncol(cellData[[feature]])]))) {
      #           outRow <- c(outRow, i)
      #         }
      #         
      #       }
      #       #remove fully nan nodes
      #       if(!is.null(outRow)) {
      #         
      #         for(feat in 1:length(cellData)) {
      #            cellData[[feat]] <- cellData[[feat]][-outRow,]
      #           
      #         }
      #       }
      
      #       #what can't be regressed put to 0 - we probably don't care anyway
      #       cellData[[feature]][is.na(cellData[[feature]])] <- 0
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
  
  
  
  #calculate mean and sd for every node
  
  for(node in (1:node_number)) {
    for(kpi in (1:kpi_number)) {
      
      data <- cellData[[kpi]][node, ]
      des_data <- cellData[[kpi]][node,]
      #seasonality
      for(hour in (0:(period-1))) {
        svec <- seq(1 + hour, length(data), period)
        festive_ind <- which(festive == T)
        
        cell_F_seasonality[[node]][kpi,hour + 1] <- mean(data[svec[festive_ind]], na.rm = T)
        cell_F_seasonality_sd[[node]][kpi,hour + 1] <- sd(data[svec[festive_ind]], na.rm = T)
        
        cell_seasonality[[node]][kpi,hour + 1] <- mean(data[svec[-festive_ind]], na.rm = T)
        cell_seasonality_sd[[node]][kpi,hour + 1] <- sd(data[svec[-festive_ind]], na.rm = T)
        
        des_data[svec[-festive_ind]] <- des_data[svec[-festive_ind]] - cell_seasonality[[node]][kpi,hour + 1] 
        des_data[svec[festive_ind]] <- des_data[svec[festive_ind]] - cell_F_seasonality[[node]][kpi,hour + 1] 
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
      
      cell_cdf[[node]][[kpi]] <- kCDF(des_data)
      
      
      #       season <- stl(times, s.window = "periodic")
      #       Fseason <- stl(Ftimes, s.window = "periodic")
      #       
      #       cell_full_season <- season$time.series[1:24]
      #       cell_F_full_season <- Fseason$time.series[1:24]
      #       
      #       cell_STL_seasonality[[node]][kpi, ] <- (cell_full_season + mean(season$time.series[,"trend"]))
      #       cell_F_STL_seasonality[[node]][kpi, ] <- (cell_F_full_season + mean(Fseason$time.series[,"trend"]))
    }
    
    
  }
  
  
  for(feature in 1:length(cellData)) {
    
    #what can't be regressed put to 0 - we probably don't care anyway
    cellData[[feature]][is.na(cellData[[feature]])] <- 0
    cellData[[feature]][is.nan(cellData[[feature]])] <- 0
    
  }  
  
  
  
  for(time in (1:final_timeslot)) {
    
    print(paste("Time: ", time))
    
    node_votes <- rep(0, node_number)
    
    
    ## START NODE OPERATION ##
    
    
    for(node in (1:node_number)) {
      
      #check time slot
      
      curr_period <- ((time-1) %% period) + 1
      
      #check if festive
      curr_day <- floor((time-1)/period) + 1
      festivity <- festive[curr_day]
      
      
      curr_values <- vector("numeric", length = kpi_number)
      
      
      if(selectMethod[1] == TRUE) {
        season_fit <- vector("numeric", length = kpi_number)
      }
      if(selectMethod[2] == TRUE) {
        distance_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[3] == TRUE) {
        trend_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[4] == TRUE) {
        previous_fit <- vector("numeric", length = kpi_number)
        
      }
      if(selectMethod[5] == TRUE) {
        history_fit <- vector("numeric", length = kpi_number)
      }
      
      
      
      #for each feature in the node, calculate measures
      for(kpi in (1:kpi_number)) {
        
        curr_values[kpi] <- cellData[[kpi]][node, time]
        
        if(festivity) {
          curr_seasonality <- cell_F_seasonality[[node]][kpi,]
          curr_seasonality_sd <- cell_F_seasonality_sd[[node]][kpi,]
        }
        else {
          curr_seasonality <- cell_seasonality[[node]][kpi,]
          curr_seasonality_sd <- cell_seasonality_sd[[node]][kpi,]
        }
        
        
        
        
        if(selectMethod[1] == TRUE) {
          season_fit[kpi] <- featureSeasonFitness(curr_seasonality,
                                                  curr_values[kpi],
                                                  curr_period, replicas = 3)
        }
        
        #find standardized values for distance and trend
        season_avg <- mean(curr_seasonality)
        season_rng <- (range(curr_seasonality)[2] 
                       - range(curr_seasonality)[1]) / 2
        
        if(is.nan(season_rng) || is.null(season_rng)) {
          print("debug")
        }
        
        if(season_rng < min_difference ) {
          season_rng <- min_difference
        }
        
        curr_val_std <- (curr_values[kpi] - season_avg) / season_rng 
        curr_season_std <- (curr_seasonality[curr_period] - season_avg) / season_rng 
        
        dist <- curr_val_std - curr_season_std
        #dist <- curr_values[kpi] - curr_seasonality[curr_period]
        
        if(curr_period > 1) {
          last_period <- curr_period - 1 
        }
        else {
          last_period <- period
        }
        
        if(time >= 2) {
          
          prev_val_std <- (cellData[[kpi]][node,time-1] - season_avg) / season_rng 
          prev_season_std <- (curr_seasonality[last_period] - season_avg) / season_rng 
          
          residual <- (curr_val_std - prev_val_std) - (curr_season_std - prev_season_std)
          #residual <- curr_values[kpi] - cellData[[kpi]][node,time-1] -(curr_seasonality[curr_period] - curr_seasonality[last_period])
          
          last_descent_mean <- cell_descent_historical_mean[[node]][kpi]
          
          cell_descent_historical_mean[[node]][kpi] <- (cell_descent_historical_mean[[node]][kpi] * (time - 1) +
                                                          residual) / time
          cell_descent_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_descent_historical_sd[[node]][kpi] ^2) ) 
                                                           + (residual - last_descent_mean) * (residual - cell_descent_historical_mean[[node]][kpi])) / time)
          
          if(is.nan(cell_descent_historical_sd[[node]][kpi])) {
            cell_descent_historical_sd[[node]][kpi] <- 0
          }
        }
        
        last_distance_mean <- cell_distance_historical_mean[[node]][kpi]
        
        cell_distance_historical_mean[[node]][kpi] <- (cell_distance_historical_mean[[node]][kpi] * (time - 1) +
                                                         dist) / time
        cell_distance_historical_sd[[node]][kpi] <- sqrt((((time - 1) * (cell_distance_historical_sd[[node]][kpi] ^2) ) 
                                                          + (dist - last_distance_mean) * (dist - cell_distance_historical_mean[[node]][kpi])) / time)
        
        
        if(is.nan(cell_distance_historical_sd[[node]][kpi])) {
          cell_distance_historical_sd[[node]][kpi] <- 0
        }
        
        if(selectMethod[4] == TRUE) {
          previous_fit[kpi] <- featureLastValuesFitness(cell_memories[[node]][kpi,],
                                                        curr_values[kpi], visual)          
        }
        
        if(cell_historical_sd[[node]][kpi] < min_difference) {
          cell_historical_sd[[node]][kpi] <- min_difference
        }
        
        if(selectMethod[2] == TRUE) {
          distance_fit[kpi] <- featureDistanceThresholdPVal(curr_season_std,
                                                            curr_val_std,
                                                            cell_distance_historical_sd[[node]][kpi],
                                                            threshold_confidence, visual)
        }
        
        if(selectMethod[3] == TRUE) {
          
          trend_fit[kpi] <- featureTrendThresholdPVal(curr_season_std,
                                                      prev_season_std,
                                                      curr_val_std,
                                                      prev_val_std,
                                                      cell_descent_historical_sd[[node]][kpi],
                                                      threshold_confidence, visual)
          
        }
        
        if(selectMethod[5] == TRUE) {
          
          history_fit[kpi] <- featureHistoryFitness(cell_historical_mean[[node]][kpi],
                                                    cell_historical_sd[[node]][kpi],
                                                    curr_values[kpi], visual)
        }
        
      }
      
    }
    
    ## END NODE OPERATION ##
    
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
                    distance_mean = cell_distance_historical_mean,
                    distance_sd = cell_distance_historical_sd,
                    kernel_cdf = cell_cdf)
  
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
featureLastValuesFitness <- function(feat_cell_past_values, feat_curr_value, visual = FALSE) {
  
  #assumes last value in feat_cell_past_values is at the end
  
  library(forecast)
  
  if(abs(max(feat_cell_past_values) - min(feat_cell_past_values)) < 0.0001) {
    #nothing to predict, all equal
    
    fcasted <- max(feat_cell_past_values)
    sdev <- 0.0001
    
  }
  else {
    
    
    #fit an arima to the historical data
    #aa <- auto.arima(feat_cell_past_values)
    
    #aa <- Arima(feat_cell_past_values, c(6,0,1))
    aa <- Arima(feat_cell_past_values, c(0,1,0), c(0,1,0), method = "CSS-ML")
    
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


featureDistanceThresholdPVal <- function(feat_cell_season_val, feat_curr_value, distance_threshold, threshold_confidence, visual = FALSE) {
  
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
  
  myDiff <- diff(season_residual, lag = lag, difference)
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


from_cdf <- function(number, cdf) {
  
  bef <- max(which(cdf$x <= number))
  aft <- min(which(cdf$x >= number))
  
  if(bef < 1) {
    p_val <- 0
  }
  else if (aft > length(cdf$x)) {
    p_val <- 1
  }
  else {
    
    p_dist <- cdf$Fhat[aft] - cdf$Fhat[bef]
    x_dist <- cdf$x[aft] - cdf$x[bef]
    
    nearB <- (number - cdf$x[bef]) / x_dist
    
    p_val <- cdf$Fhat[bef] + p_dist * (nearB)
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

fastGeneralWeightedFisherScan <- function(priorities, weights, anomThreshold = 0.95, base_priority = 0.68) {
  
  #transform priorities according to their weight
  for(p in (1:length(priorities))) {
    
    priorities[p] <- priorityWeight(priorities[p], weights[p], base_priority)
    
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
  
  #Best result
  
  bestSize <- 0
  bestScore <- 0
  
  #do not report clusters of less than 3 nodes
  minSize <- 3
  
  for(size in c(minSize:length(prior_order))) {
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
fisher.baseChange <- function(normPriorities) {
  
  changingb <- c(0:1000)/1000
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

fisher.tryConvergence <- function(priorities = c(0.5, 0.3, 0.7, 0.3), weights = c(0.7, 0.8, 0.8, 0.3), 
                                  baseStart = 0.7, iterations = 10) {
  
  
  baselines <- rep(baseStart, length(priorities))
  
  errorGrowth <- NULL
  histErrors <- NULL
  
  string <- paste(baselines, collapse = ", ")
  print(paste("Start baselines: ",string))
  
  
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
  
  histErrors <- c(histErrors, sum(lasterrors ^ 2))
  
  #Iterative process
  for(it in 1:iterations) {
    
    print(paste("Iteration ", it))
    
    newBaselines <- baselines
    weightPriorities <- priorities
    
    print("Computing weighted Priorities")
    for(i in 1:length(priorities)) {
      #compute weightede priority with old weights
      weightPriorities[i] <- priorityWeight(priorities[i], weights[i], baselines[i])
      
    }
    
    errors <- NULL
    print("Computing new baselines")
    for(i in 1:length(priorities)) {
      #compute baselines that minimizes difference
      fisher_test <- fisher.baseChange(weightPriorities[-i])
      newBaselines[i] <- fisher_test$base 
      errors <- c(errors, fisher_test$error)
    }
    
    errorGrowth <- errors - lasterrors #c(errorGrowth, sum(errors ^ 2))
    
    newDiff <- baselines - newBaselines
    
    growing <- which(errorGrowth > 0)
    
    newBaselines[growing] <- baselines[growing] 
    
    
    string <- paste(newBaselines, collapse = ", ")
    print(paste("New baselines: ",string))
    baselines <- newBaselines
    lasterrors <- errors
    histErrors <- c(histErrors, sum(errors ^ 2))
  }
  
  plot(0:iterations, histErrors, type = "o", main ="SSError / iteration")
  
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

#surpassed by training ?
updateSeasonality <- function(node_seasonality, curr_values, curr_period, curr_status, noSeasonPresent = FALSE) {
  
  cellUpdateRateOK <- 0.7
  cellUpdateRateOUT <- 0.3
  
  cellUpdateRate <- 0
  
  if(noSeasonPresent == TRUE) {
    
    node_seasonality[,curr_period] <- curr_values
    
  }
  else {
    
    if(curr_status == "OK") {
      
      cellUpdateRate <- cellUpdateRateOK
      
    }
    else if(curr_status == "Anomaly") {
      
      cellUpdateRate <- cellUpdateRateOUT
      
    }
    else {
      #default
      cellUpdateRate <- cellUpdateRateOK
    }
    
    node_seasonality[,curr_period] <- ((1 - cellUpdateRate) * node_seasonality[,curr_period]) + (cellUpdateRate * curr_values)
    
  }
  
  return(node_seasonality)
  
  
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

shiftRightMatrix <- function(mat, num = 1) {
  
  m1 <- matrix(mat[,(ncol(mat) - num + 1):ncol(mat)], ncol = num)
  m2 <- matrix(mat[,1:(ncol(mat)-num)], ncol = ncol(mat) - num)
  
  mat <- cbind(m1,m2)
  return(mat)
  
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



myTimeImage <- function(coord, vals, min = NULL, max = NULL,title = NULL) {
  
  if(is.null(min)){
    min <- min(vals)
  }
  if(is.null(max)) {
    max <- max(vals)
  }
  
  vals <- rowMeans(as.matrix(vals), na.rm = T)
  
  mypal <- colorRampPalette( c( "blue", "red" ) )( 5 )
  
  cols <- map2color(vals, mypal, limits = c(min, max))
  
  plot(coord[,1], coord[,2], col = cols, main = title, pch = 19)
  
  
}

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
      
      distance_matrices <- c(distance_matrices, list(dist(anomalies_by_length[[i]], method = "DTW")))
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








###################### Motif Forecast #########


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

motif_search <- function(current_data, level_set, motif_set, k) {
  
  
  lab_set <- as.matrix(label_with_interval(current_data, level_set))
  
  an_support <- 0
  
  for(motif in motif_set) {
    #Compare beginning of motif with current data
    big_dist_threshold <- k
    
    main_dist_threshold <- big_dist_threshold * 0.1
    neigh_dist_threshold <- big_dist_threshold * 0.1 * (nrow(lab_set) - 1) 
    
    candidate_motif <- as.matrix(motif$chain[,1])
    
    main_distance <- length(which(candidate_motif[1] != lab_set[1]))
    neigh_distance <- length(which(candidate_motif[-1] != lab_set[-1]))
    
    main_distance <- sum(abs(candidate_motif[1] - lab_set[1]))
    neigh_distance <- sum(abs(candidate_motif[-1] - lab_set[-1]))
    
    match_level <- 0
    
    if(main_distance <= main_dist_threshold) {
      if(neigh_distance <= neigh_dist_threshold) {
        match_level <- 1
      }
      else {
        match_level <- 0.8
      }
    }
    else if(neigh_distance <= neigh_dist_threshold) {
      match_level <- 0.4
      
    }    
    
    #if matching, calculate forecast as support of this motif, then forecast for next of this motif
    an_support <- 1 - ((1 - an_support) * (1 - match_level*motif$support))
    
  }
  
  
  return(an_support)
  
}

support_calc <- function(lab_set, chain, k) {
  
  #calculate a support measure for current candidate_motif
  #how many times this pattern was followed by correct prediction
  
  calced_chains <- NULL
  
  candidate_motif <- as.matrix(chain)
  
  window_size <- ncol(candidate_motif)
  
  big_dist_threshold <- k
  
  main_dist_threshold <- big_dist_threshold * 0.1
  neigh_dist_threshold <- big_dist_threshold * 0.2 * (nrow(lab_set) - 1) 
  
  
  support <- NULL
  
  match_times <- 0
  
  next_match <- 0
  next_support <- 0
  
  for(t in c(1:(ncol(lab_set)))) {
    
    matching <- 0
    
    #How many categorical labels are different?
    main_distance <- length(which(candidate_motif[1,1] != lab_set[1,t]))
    neigh_distance <- length(which(candidate_motif[-1,1] != lab_set[-1,t]))
    
    main_distance <- sum(abs(candidate_motif[1,1] - lab_set[1,t]))
    neigh_distance <- sum(abs(candidate_motif[-1,1] - lab_set[-1,t]))
    
    if(main_distance <= main_dist_threshold) {
      
      if(neigh_distance <= neigh_dist_threshold) {
        
        matching <- 1
        match_times <- match_times + 1
        
        
      }
      else {
        
        #matches own history, not neighborhood        
        matching <- 0.8
        match_times <- match_times + 0.8
        
      }
    }
    else if(neigh_distance <= neigh_dist_threshold) {
      
      matching <- 0.4
      match_times <- match_times + 0.4
      
    }
    
    if(matching > 0 && ncol(candidate_motif) > 1 && t < ncol(lab_set) ) {
      
      #check if next is matching
      main_distance <- length(which(candidate_motif[1,2] != lab_set[1,t+1]))
      neigh_distance <- length(which(candidate_motif[-1,2] != lab_set[-1,t+1]))
      
      main_distance <- sum(abs(candidate_motif[1,2] - lab_set[1,t+1]))
      neigh_distance <- sum(abs(candidate_motif[-1,2] - lab_set[-1,t+1]))
      
      if(main_distance <= main_dist_threshold) {
        
        if(neigh_distance <= neigh_dist_threshold) {
          
          next_match <- next_match + 1
          
          
        }
        else {
          
          #matches own history, not neighborhood      
          next_match <- next_match + 0.8
          
        }
      }
      else if(neigh_distance <= neigh_dist_threshold) {
        
        next_match <- next_match + 0.4
        
      }
    }
    
  }
  
  support <- 1
  
  #calculate next_support
  if(ncol(candidate_motif) > 1) {
    
    next_calc <- support_calc(lab_set, candidate_motif[,-1], k)
    calced_chains <- c(calced_chains, next_calc$chains)
    
    support <- c(support, (next_match / match_times) * next_calc$lastSupport)
    
  }
  
  calced_chains <- c(list(list(chain = candidate_motif,
                               support = support)), calced_chains)
  
  return(list(chains = calced_chains,
              lastSupport = support))
  
  
}

anomaly_motifs <- function (training_set, anomaly_set, min_k = 1, max_k = 24) {
  
  print("Labeling set...")
  #divide training_set in a group of labels
  labeled_set <- psf_label_set(training_set, min_k, max_k)
  
  motif_set <- NULL
  
  lab_set <- labeled_set$labels
  liv_set <- labeled_set$levels
  
  current_motif <- NULL
  current_support <- 0
  
  print("Isolating meaningful motifs...")
  
  for(cand in (anomaly_set)) {
    
    #find motif representing found anomaly
    current_motif <- cbind(current_motif, 
                           as.matrix(lab_set[, (cand)]))
    
    if( !((cand + 1) %in% anomaly_set)) {
      motif_set <- c(motif_set, support_calc(lab_set, current_motif, max_k)$chains)
      current_motif <- NULL
      
    }
  }
  
  
  return(list(set = motif_set,
              lev = liv_set))
  
}

motif_anomaly_discovery <- function (training_set, anomaly_set, n, offset = 1, min_window = 1, max_window = 24, min_k = 1, max_k = 24) {
  
  
  
  motif_set <- NULL
  candidate_set <- NULL
  outlier_set <- NULL
  support_set <- NULL
  
  result <- NULL
  
  
  print("Labeling set...")
  #divide training_set in a group of labels
  labeled_set <- psf_label_set(training_set, min_k, max_k)
  
  lab_set <- labeled_set$labels
  liv_set <- labeled_set$levels
  
  print("Finding best window...")
  #find best window to use for forecasting
  window_disc <- bestWindowDiscovery(training_set, liv_set, n, offset, min_window, max_window, min_k, max_k, return_error = TRUE)
  window_size <- window_disc$win_size
  
  
  big_dist_threshold <- max_k * window_size
  
  main_dist_threshold <- big_dist_threshold * 0.1
  neigh_dist_threshold <- big_dist_threshold * 0.2 * (nrow(training_set) - 1) 
  
  
  
  print("Isolating meaningful motifs...")
  for(cand in (anomaly_set)) {
    
    #find motifs (windows) before candidate_set.
    ##find windows of nn too
    candidate_motif <- as.matrix(lab_set[, (cand - offset - window_size +1):(cand - offset)])
    
    
    matches <- 0
    confirmed <- TRUE
    
    for(t in c((window_size + offset):ncol(lab_set))) {
      
      #How many categorical labels are different?
      main_distance <- length(which(candidate_motif[1,] != lab_set[1,((t - window_size - offset + 1) : (t - offset))]))
      neigh_distance <- length(which(candidate_motif[-1,] != lab_set[-1,((t - window_size - offset + 1) : (t - offset))]))
      
      main_distance <- sum(abs(candidate_motif[1,] - lab_set[1,((t - window_size - offset + 1) : (t - offset))]))
      neigh_distance <- sum(abs(candidate_motif[-1,] - lab_set[-1,((t - window_size - offset + 1) : (t - offset))]))
      
      if(main_distance <= main_dist_threshold) {
        
        #found match
        #If the motif only appear before candidates, consider it a valid motif (with given support), if not discard it
        ##check motif of nn too
        
        
        if(neigh_distance <= neigh_dist_threshold) {
          
          if(!(t %in% anomaly_set)) {
            #motif appears before normal data
            #confirmed <- FALSE
            matches <- matches - 1
            
          }
          else {
            #perfect match on another candidate
            #increase support
            matches <- matches + 1
          }
          
        }
        else {
          
          #matches own history, not neighborhood        
          
          if(!(t %in% anomaly_set)) {
            #motif appears before normal data
            #confirmed <- FALSE
            matches <- matches - 0.8
            
          }
          else {
            #close match on another candidate
            #increase support
            matches <- matches + 0.8
          }
          
        }
        
      }
      else if(neigh_distance <= neigh_dist_threshold) {
        
        
        #matches neigh history, not own... outlier?            
        if(!(t %in% anomaly_set)) {
          #motif appears before normal data
          matches <- matches - 0.4
          
        }
        else {
          #far match on another candidate
          #increase support
          matches <- matches + 0.4
        }
        
      }
    }
    
    if(confirmed == TRUE) {
      
      outlier_set <- c(outlier_set, cand)
      motif_set <- c(motif_set, list(candidate_motif))
      support_set <- c(support_set, matches)
      
      result <- c(result, list(motif = list(candidate_motif),
                               support = matches))
    }
  }
  
  
  return(result)
  
  
  
}


###################### Linked List Forecast #############################


compareHistory <- function(vec1, vec2) {
  
  
  vecDiff <- vector("numeric", length = length(vec1))
  
  #vecDiff <- abs(vec1 - vec2)
  
  
  for(i in (1:length(vec1))) {
    
    if(is.na(vec1[[i]]) && is.na(vec2[[i]])) {
      #where v1 and v2 NA, consider distance 1
      vecDiff[i] <- 0.5
    }
    else if(is.na(vec1[[i]]) || is.na(vec2[[i]])) {
      #where only one of v1 and v2 is NA, consider distance 2 
      vecDiff[i] <- 0.5
    }
    else {
      vecDiff[i] <- compareStates(vec1[[i]], vec2[[i]])
    }
  }
  
  
  #where v1 and v2 NA, consider distance 1
  #vecDiff[is.na(vec1) & is.na(vec2)] <- 1
  
  #where only one of v1 and v2 is NA, consider distance 2
  #vecDiff[is.na(vec1) != is.na(vec2)] <- 3
  
  comparison <- max(0, mean(vecDiff))
  
  return(comparison)
  
}


testCompareHistory <- function(vec1, vec2) {
  
  
  
  vecDiff <- abs(vec1 - vec2)
  
  
  
  #where v1 and v2 NA, consider distance 1
  vecDiff[is.na(vec1) & is.na(vec2)] <- 1
  
  #where only one of v1 and v2 is NA, consider distance 2
  vecDiff[is.na(vec1) != is.na(vec2)] <- 3
  
  comparison <- max(0, 1 - mean(vecDiff^2) / 10)
  
  return(comparison)
  
}

#compares two anomaly_table[[node]][[anomaly]]$precedent or anomaly_table[[node]][[anomaly]]$subsequent
compareStates <- function(list1, list2) {
  
  
  values1 <- vector("numeric", length = length(list1))
  values2 <- vector("numeric", length = length(list2))
  
  fac1 <- vector("character",length = length(list1))
  fac2 <- vector("character",length = length(list2))
  
  for(i in (1:length(list1))) {
    
    values1[i] <- list1[[i]]$value
    values2[i] <- list2[[i]]$value
    
    
    
    fac1[i] <- list1[[i]]$state
    fac2[i] <- list2[[i]]$state
    
  }
  
  
  ord1 <- as.numeric(ordered(fac1, levels = c("OUT-down", "TREND-down", "OK", "TREND-up", "OUT-up")))
  ord2 <- as.numeric(ordered(fac2, levels = c("OUT-down", "TREND-down", "OK", "TREND-up", "OUT-up")))
  
  statDist <- abs(ord2 - ord1)
  #transform statDist from range [0,4] to range [0.5, 1.5]
  statDist <- (statDist / 4) + 0.5
  
  #valDist <- ((abs(values2 - values1) + 1)^2) - 1
  valDist <- (abs(values2 - values1))
  
  #weightedDist <- statDist * valDist
  weightedDist <- valDist
  
  comparison <- max(0, 1 - sqrt(mean(weightedDist)))
  
  return(comparison)
  
  
}



#input: anomaly structure ( for a single node ), the [look_behind] preceding elements, the current element
#returns: the probability of reaching an anomalous state within look_ahead prediction
predictAnomaly <- function(anomalyStructure, curr_state, prec_states, look_ahead = 3, genericPrediction = FALSE) {
  
  curr_hops <- NULL
  next_hops <- NULL
  anom_hops <- NULL
  
  
  #initialize start of search
  curr_entry_id <- structure.searchAnomaly(anomalyStructure, curr_state)
  
  if(!is.null(curr_entry_id)) {
    
    curr_entry <- structure.getEntryById(anomalyStructure, curr_entry_id)
    
    #hopProb = probability of being in that hop \ 1 for current state
    curr_hops <- c(curr_hops, list(list(hopState = curr_entry,
                                        hopProb = 1)))
    
    for(step in c(0:look_ahead)) {
      
      next_hops <- NULL
      
      for(hop in curr_hops) {
        
        if(hop$hopState$flag == "ANOMALY") {
          #add anomaly reachable in [step] time, with probability hop$hopProb
          anom_hops <- c(anom_hops, list(c(hop, steps = step)))
        }
        else{
          
          identifiedLine <- NULL
          #current similarity measure
          similarHist <- 0
          
          #do not consider at all lines less similar than this
          similarThreshold <- 0.5
          
          #search in the state the next elements corresponding to precedingIds
          for(line in hop$hopState$next_mat) {
            
            
            
            similarity <- compareHistory(line$history, prec_states)
            #choose line whose history match best current one          
            if(similarity > max(similarHist, similarThreshold)) {
              #found line better than current one
              
              similarHist <- similarity
              identifiedLine <- line
              
            } 
            
          }
          #add next hops to list
          if(!is.null(identifiedLine)) {
            for(nhop in identifiedLine$jumps) {
              
              probb <- hop$hopProb * nhop$prob * similarHist
              next_hops <- c(next_hops, list(list(hopState = getState(anomalyStructure, nhop$nextState),
                                                  hopProb = (hop$hopProb * nhop$prob * similarHist))))
              
              wait <- 1
            }
          }
          
        }
        
      }
      
      curr_hops <- next_hops
      
      if(genericPrediction == T) {
        
        for(nhop in next_hops) {
          if(nhop$hopProb > 0.8)
            print(paste("Next prediction ",nhop$hopState$id, ", prob ",nhop$hopProb))
          
        }
      }
      
    }
    
  }
  #anom_hops contains all reached anomalies
  
  anomProbability <- 0
  for(anomaly in anom_hops) {
    
    anomProbability <- anomProbability + anomaly$hopProb
    
  }
  
  return(anom_hops)
  
  
  
  
}

testPredictAnomaly <- function(anomStruc) {
  
  
  look_behind <- 4
  range <- 15
  precState <- c(NA,NA,NA,NA)
  
  #generate initial predecessors
  for(t in (1:look_behind)) {
    
    mSnapshot <- runif(1, 1,range)
    stateId <- round(mSnapshot)
    
    precState[t] <- stateId
    
    
  }
  
  for(i in c(1:500)) {
    
    mSnapshot <- runif(1, 1,range)
    stateId <- round(mSnapshot)
    
    predAnom <- predictAnomaly(anomStruc, stateId, precState,look_ahead = 2)
    
    if(predAnom >0) {
      
      print(paste("stateId ",stateId, "anomalous w story ",paste(precState[], collapse = " "), " with prob ",predAnom))
      
    }
    
    
    precState <- shiftRight(precState)
    precState[1] <- stateId
    
  }
  
  
  
}



createTestAnomalyStruc <- function() {
  
  anomalyStruct <- NULL
  
  #precState[1] last hop, precState[2] second to last etc...
  precState <- c(NA,NA,NA,NA)
  
  look_behind <- 4
  initialProb <- 0.8
  range <- 15
  
  #generate ids of anomalous states
  anomalyIds <- sample(1:range, 3)
  
  #generate initial predecessors
  for(t in (1:look_behind)) {
    
    mSnapshot <- runif(1, 1,range)
    stateId <- round(mSnapshot)
    #regenerate if anomaly
    while(any(stateId == anomalyIds[], na.rm=T)) {
      mSnapshot <- runif(1, 1,range)
      stateId <- round(mSnapshot)
    }
    
    precState[t] <- stateId
    
    
    
  }
  
  
  for(i in c(1:250)) {
    
    #snapshot(state val) -> stateId
    
    mSnapshot <- runif(1, 1,range)
    
    stateId <- round(mSnapshot)
    
    
    if(any(stateId == anomalyIds[], na.rm=T)) {
      #add anomaly if not present
      
      if(is.null(getState(anomalyStruct, stateId))) {
        anomalyState <- list(list(id = stateId,
                                  snapshot = mSnapshot,
                                  next_mat = NULL,
                                  flag = "ANOMALY"))
        
        anomalyStruct <- c(anomalyStruct, anomalyState)
      }
      
      #add predecessors
      
      curr_st <- stateId
      
      for(ind in (1:length(precState))) {
        
        predecessor <- precState[ind]
        beforePredecessor <- NULL
        
        if(ind < length(precState)) {
          beforePredecessor <- c(precState[(ind + 1):length(precState)], 
                                 rep(NA, ind))
        }
        else{
          beforePredecessor <- c(rep(NA, ind))
        }
        
        #number of valid states memorized before predecessor
        beforeCount <- sum(!is.na(beforePredecessor))
        
        predecessorEntry <- getState(anomalyStruct, precState[ind])
        
        
        if(is.null(predecessorEntry)) {
          
          #if predecessor is not in structure
          predecessorEntry <- list(list(id = predecessor,
                                        snapshot = predecessor,
                                        next_mat = list(list(history = beforePredecessor,
                                                             jumps = list(list(nextState = curr_st,
                                                                               prob = initialProb)))),
                                        flag = "TREND"))
          
          
          anomalyStruct <- c(anomalyStruct, predecessorEntry)
        }
        else {
          
          #if predecessor is in structure
          #find index of predecessor within structure
          predecessorIndex <- getStateIndex(anomalyStruct, precState[ind])
          
          #check if have to add new line
          jum <- getJumps(predecessorEntry, beforePredecessor)
          
          if(is.null(jum)) {
            #add new line
            line <- list(list(history = beforePredecessor,
                              jumps = list(list(nextState = curr_st,
                                                prob = initialProb)))) 
            
            predecessorEntry$next_mat <- c(predecessorEntry$next_mat,
                                           line)
            
          }
          else {
            
            #already present jumps for same history
            predecessorLineIndex <- getJumpsIndex(predecessorEntry, beforePredecessor)
            
            #check if curr_st is in nexts for the line
            
            present <- F
            for(nextInd in (1:length(predecessorEntry$next_mat[[predecessorLineIndex]]$jumps))) {
              if(curr_st == jum[[nextInd]]$nextState) {
                #increment probability
                jum[[nextInd]]$prob <- min(1, jum[[nextInd]]$prob + 0.3)
                present = T 
              }
              else {
                #decrement probability
                jum[[nextInd]]$prob <- max(0, jum[[nextInd]]$prob - 0.2)
                
              }
            }
            
            #if not, add state to line
            if(!present) {
              
              state <- list(list(nextState = curr_st,
                                 prob = 0.8/length(jum)))
              
              jum <- c(jum, state)
              
            }
            
            #update state
            predecessorEntry$next_mat[[predecessorLineIndex]]$jumps <- jum
            
            
          }
          
          #update struct 
          anomalyStruct[[predecessorIndex]] <- predecessorEntry
          
          
        }
        
        #proceed with the chain
        curr_st <- precState[ind]
        
      }
      
    }
    else {
      
      #not anomaly, check in states for matching history
      
      st <- getState(anomalyStruct,stateId)
      if(!is.null(st)) {
        
        
        
      }
      
      
    }
    
    precState <- shiftRight(precState)
    precState[1] <- stateId
    
  }
  
  
  
  return(anomalyStruct)
  
  
  
}




#input: anomaly structure for a node
#returns: the id state, or NULL if not present
getState <- function(anomalyStructure, id) {
  
  for(state in anomalyStructure) {
    
    if(state$id == id)
      return(state)
  }
  
  return(NULL)
  
}

getStateIndex <- function(anomalyStructure, id) {
  
  for(j in (1:length(anomalyStructure))) {
    
    state <- anomalyStructure[[j]]
    
    if(state$id == id)
      return(j)
  }
  
  return(NULL)
  
}

getJumps <- function(anomalyEntry, predecessorsIds) {
  
  
  #search in the state the next elements corresponding to precedingIds
  for(line in anomalyEntry$next_mat) {
    
    if(identical(line$history, predecessorsIds)) {
      #found line
      
      return(line$jumps)
      
    } 
  }
  
  return(NULL)
}

getJumpsIndex <- function(anomalyState, predecessorsIds) {
  
  #search in the state the next elements corresponding to precedingIds
  for(j in (1:length(anomalyState$next_mat))) {
    
    line <- anomalyState$next_mat[[j]]
    
    if(identical(line$history, predecessorsIds)) {
      #found line
      
      return(j)
      
    } 
  }
  
  return(NULL)
}

getNexts <- function(anomalyState, predecessorsIds) {
  
  #search in the state the next elements corresponding to precedingIds
  for(line in anomalyState$next_mat) {
    
    if(identical(line$history,predecessorsIds)) {
      #found line
      
      return(line$jumps)
      
    } 
  }
  
  return(NULL)
}


#returns the id of an entry with most similar state, or NULL
#if sinks is TRUE may return entries even if not anomalies or trending to anomalies
structure.searchAnomaly <- function(anomalyStruct, state, sinks = FALSE) {
  
  #current similarity measure
  bestSimilarity <- 0
  
  similarStateThreshold <- 0.8
  
  identifiedId <- NULL
  
  #do not consider at all entries less similar than this
  
  for(entry in anomalyStruct) {
    
    if(sinks == TRUE || entry$flag != "SINK") {
      
      similarity <- compareStates(state, entry$snapshot)
      
      
      #choose best entry         
      if(similarity > max(bestSimilarity, similarStateThreshold)) {
        #found line better than current one
        
        bestSimilarity <- similarity
        identifiedId <- entry$id
        
      }
      
    }
    
  }
  
  return(identifiedId)
  
}

structure.getEntryById <- function(anomalyStruct, id) {
  
  
  for(j in (1:length(anomalyStruct))) {
    
    entry <- anomalyStruct[[j]]
    
    if(entry$id == id)
      return(entry)
  }
  
  return(NULL)
  
}

structure.getIndexById <- function(anomalyStruct, id) {
  
  
  for(j in (1:length(anomalyStruct))) {
    
    state <- anomalyStruct[[j]]
    
    if(state$id == id)
      return(j)
  }
  
  return(NULL)
  
}

#returns an unused id
structure.newId <- function(anomalyStruct) {
  
  anomaly_max_id <- 1
  
  for(j in (1:length(anomalyStruct))) {
    
    anomaly_max_id <- max(anomaly_max_id, anomalyStruct[[j]]$id)
    
  }
  
  anomaly_unusedId <- anomaly_max_id + 1
  
  return(anomaly_unusedId)
  
}

#returns the INDEX of the line in an entry with best similar history, or nothing
entry.searchLine <- function(entry, prec_history) {
  
  #current similarity measure
  bestSimilarity <- 0
  
  similarHistoryThreshold <- 0.7
  
  identifiedLine <- NULL
  identifiedInd <- NULL
  if(!is.null(entry$next_mat)) {
    #search in the state the next elements corresponding to precedingIds
    for(i in (1:length(entry$next_mat))) {
      
      line <- entry$next_mat[[i]]
      
      similarity <- compareHistory(line$history, prec_history)
      
      #choose line whose history match best current one          
      if(similarity > max(bestSimilarity, similarHistoryThreshold)) {
        #found line better than current one
        
        bestSimilarity <- similarity
        identifiedLine <- line
        identifiedInd <- i
        
      } 
      
    }
  }
  
  
  return(identifiedInd)
  
}

#add a new anomaly with curr_state and prec_states
#returns an updated copy of the structure
structure.addNewAnomaly <- function(anomalyStruct, curr_state, prev_states) {
  
  curr_state_id <- structure.searchAnomaly(anomalyStruct, curr_state)
  
  
  #if not present, add new entry for anomaly state
  if(is.null(curr_state_id)) {
    
    curr_state_id <- structure.newId(anomalyStruct)
    
    anomalyEntry <- list(id = curr_state_id,
                         snapshot = curr_state,
                         next_mat = NULL,
                         flag = "ANOMALY")
    
    anomalyStruct <- c(anomalyStruct, list(anomalyEntry))
  }
  
  
  print(paste("Current state ",curr_state_id))
  
  #add predecessors
  
  for(ind in (1:length(prev_states))) {
    
    predecessor <- prev_states[[ind]]
    beforePredecessor <- NULL
    
    if(ind < length(prev_states)) {
      beforePredecessor <- c(prev_states[(ind + 1):length(prev_states)], 
                             rep(NA, ind))
    }
    else{
      beforePredecessor <- c(rep(NA, ind))
    }
    
    
    anomalyStruct <- structure.updateEntry(anomalyStruct, curr_state_id, predecessor, 
                                           beforePredecessor, addNewEntry = TRUE)
    
    predecessorId <- structure.searchAnomaly(anomalyStruct, predecessor, sinks = TRUE)                   
    
    #proceed with the chain
    curr_state_id <- predecessorId
    
    
    print(paste("Current state ",curr_state_id))
    
  }
  
  
  return(anomalyStruct)
  
}

structure.normalUpdate <- function(anomalyStruct, curr_state, prev_states) {
  
  ind <- 1
  
  predecessor <- prev_states[[ind]]
  beforePredecessor <- NULL
  
  if(ind < length(prev_states)) {
    beforePredecessor <- c(prev_states[(ind + 1):length(prev_states)], 
                           rep(NA, ind))
  }
  else{
    beforePredecessor <- c(rep(NA, ind))
  }
  
  predecessorId <- structure.searchAnomaly(anomalyStruct, predecessor, sinks = FALSE)
  
  
  
  if(!is.null(predecessorId)) {
    
    #search sink state
    curr_state_id <- structure.searchAnomaly(anomalyStruct, curr_state, sinks = T)
    
    #if not present, add new entry for sink state
    if(is.null(curr_state_id)) {
      
      curr_state_id <- structure.newId(anomalyStruct)
      
      sinkEntry <- list(id = curr_state_id,
                        snapshot = curr_state,
                        next_mat = NULL,
                        flag = "SINK")
      
      anomalyStruct <- c(anomalyStruct, list(sinkEntry))
      
    }
    
    
    print(paste("Current state ",curr_state_id))
    
    #update entry of predecessor
    anomalyStruct <- structure.updateEntry(anomalyStruct, curr_state_id, predecessor, 
                                           beforePredecessor, addNewEntry = FALSE)
    
  }
  
  
  
  return(anomalyStruct)
}

#addNew means that the current update is done when an anomaly is detected (and lines and entries must be added)
structure.updateEntry <- function(anomalyStruct, curr_state_id, predecessor, beforePredecessor, addNewEntry = TRUE) {
  
  #search if previous was trending to anomaly, if addNew = FALSE
  predecessorId <- structure.searchAnomaly(anomalyStruct, predecessor, sinks = addNewEntry)
  
  if(is.null(predecessorId) && addNewEntry == TRUE) {
    
    predecessorId <- structure.newId(anomalyStruct)
    
    #if predecessor is not in structure
    predecessorEntry <- list(id = predecessorId,
                             snapshot = predecessor,
                             next_mat = NULL,
                             flag = "TREND")
    
    predecessorEntry <- entry.updateLine(predecessorEntry, curr_state_id, 
                                         beforePredecessor, addNewLine = TRUE)
    
    
    anomalyStruct <- c(anomalyStruct, list(predecessorEntry))
    
    
  }
  else if(!is.null(predecessorId)) {
    
    #if predecessor is in structure as trending
    #find index of predecessor within structure
    predecessorIndex <- structure.getIndexById(anomalyStruct, predecessorId)
    predecessorEntry <- anomalyStruct[[predecessorIndex]]
    
    #updateEntry with new lines
    predecessorEntry <- entry.updateLine(predecessorEntry, curr_state_id, 
                                         beforePredecessor, addNewLine = addNewEntry)
    
    if(predecessorEntry$flag == "SINK")
      predecessorEntry$flag = "TREND"
    
    anomalyStruct[[predecessorIndex]] <- predecessorEntry
  }
  
  return(anomalyStruct)
  
}

#returns updated predecessorEntry
entry.updateLine <- function(predecessorEntry, curr_state_id, beforePredecessor, addNewLine = TRUE) {
  
  initialProb <- 0.8
  initialValidity <- 1
  
  updateCorrectRate <- 0.5
  updateIncorrectRate <- 0.3
  
  similarHistoryThreshold <- 0.7
  
  
  #check if have to add new line
  predecessorLineIndex <- entry.searchLine(predecessorEntry, beforePredecessor)
  
  if(is.null(predecessorLineIndex) && addNewLine == TRUE) {
    
    #add new line
    
    line <- list(history = beforePredecessor,
                 jumps = list(list(nextState = curr_state_id,
                                   prob = initialProb)),
                 validity = initialValidity)
    
    predecessorEntry$next_mat <- c(predecessorEntry$next_mat,
                                   list(line))
    
  }
  else if(!is.null(predecessorLineIndex)) {
    
    #already present line or lines for same history
    
    for(lineIndex in  c(1:length(predecessorEntry$next_mat))) {
      
      line <- predecessorEntry$next_mat[[lineIndex]]
      similarity <- compareHistory(line$history, beforePredecessor)
      
      #only update lines somewhat similar
      if(similarity > similarHistoryThreshold) {
        
        #found similar line
        #start updating
        
        line$validity <- line$validity + 1
        jum <- line$jumps
        
        #check if curr_st is in nexts for the line
        
        present <- F
        totalProb <- 0
        
        for(nextInd in (1:length(jum))) {
          if(curr_state_id == jum[[nextInd]]$nextState) {
            #increment probability
            jum[[nextInd]]$prob <- jum[[nextInd]]$prob + similarity*updateCorrectRate
            present = T 
          }
          else {
            #decrement probability
            jum[[nextInd]]$prob <- max(0,jum[[nextInd]]$prob - similarity*updateIncorrectRate)
            
          }
          totalProb <- totalProb + jum[[nextInd]]$prob
        }
        
        #if not, add nextstate to line
        if(!present) {
          
          state <- list(nextState = curr_state_id,
                        prob = (similarity*initialProb)) #/line$validity
          
          totalProb <- totalProb + state$prob
          
          jum <- c(jum, list(state))
          
        }
        
        
        #set probability to sum = 1
        for(nextInd in (1:length(jum))) {
          if(totalProb != 0) {
            jum[[nextInd]]$prob <- jum[[nextInd]]$prob / totalProb
          }
        }
        
        line$jumps <- jum
        
        predecessorEntry$next_mat[[lineIndex]] <- line
        
      } 
      
    }
    
  }
  
  return(predecessorEntry)
  
  
}


buildStateEntry <- function(node, prev_state, prev_neigh_value, state_mem, cluster_size) {
  
  entry <- NULL
  
  #entry[1:t]
  for(t in (1:state_mem)) {
    
    this_l <- NULL
    
    #this_l[[time]][nodes]
    for(i in (1:cluster_size)) {
      this_l <- c(this_l, list(list(state = prev_state[[node]][i,t], 
                                    value = prev_neigh_value[[node]][i,t])))
    }
    
    
    entry <- c(entry, list(this_l))
  }
  
  return(entry)
  
  
}

buildStateEntrySeason  <- function(node, prev_state, prev_neigh_value, state_mem, cluster_size, voteResult, weights) {
  
  entry <- NULL
  
  #entry[1:t]
  for(t in (1:state_mem)) {
    
    this_l <- NULL
    val <- 0
    #this_l[[time]][nodes]
    
    
    val <- mean(prev_neigh_value[[node]][,t] * weights)
    
    this_l <- c(this_l, list(list(state = voteResult, 
                                  value = val)))
    
    entry <- c(entry, list(this_l))
  }
  
  return(entry)
  
  
}


###################### Markov Previsions ###################################

#returns a list containing next step_numb steps with attached probability according to the cell historical data given
#MST = Markov State Transition (extended with eventual history)
#cellFeatures: a list of n arrays, each containing features value in a window of time
cellPrevisionMarkov <- function(cellFeatures, cellMST, step_numb) {
  
  predictedFeatures <- rep(NULL, length(cellFeatures))
  
  for(featureIndex in c(1:length(cellFeatures))) {
    
    featureValues <- cellFeatures[[featureIndex]]
    cellFeatMST <- cellMST[[featureIndex]]
    
    predictedFeatures[[featureIndex]] <- c(predictedFeatures, 
                                           featurePrevisionMarkov(featureValues, cellFeatMST, step_numb))
    
    
    
  }
  
  return(predictedFeatures)
  
}

maxiTestUpdate <- function() {
  
  feat_numb <- 1
  memory <- 5
  CSMT <- rep(list(NULL), feat_numb)
  prev_states <- rep(list(rep(NA, memory)), feat_numb)
  
  testData <- NULL
  predData <- NULL
  
  for(i in c(1:500)) {
    
    #data reception etc...
    for(j in c(1:feat_numb)) {
      
      rand_state <- round(i %% 5 + runif(1,  -0.7,+0.7))
      
      testData <- c(testData, rand_state)
      
      feature_hist <- prev_states[[j]]
      #update history
      feature_hist <- c(rand_state, feature_hist[-length(feature_hist)])
      prev_states[[j]] <- feature_hist
      
      
    }
    
    #data prediction (just update here)
    CSMT <- updateCellMarkov(prev_states, CSMT)
    
    predictedVals <- cellPrevisionMarkov(prev_states, CSMT, 2)
    
    bestPrediction <- NULL
    bestProb <- 0
    
    for(prediction in predictedVals[[1]]) {
      
      if(prediction$prob > bestProb) {
        bestPrediction <- prediction
        bestProb <- prediction$prob
      }
      
    }
    
    if(is.null(bestPrediction)) {
      predData <- c(predData, NA)
    }
    else{
      predData <- c(predData, bestPrediction$curr_val)
    }
  }
  
  plot(testData[450:500])
  lines(testData[450:500], col = "blue")
  points(predData[450:500], col="red")
  lines(predData[450:500], col = "red")
  
  return(CSMT)
  
  
}



updateCellMarkov <- function(cellFeatures, cellMST) {
  
  for(featureIndex in c(1:length(cellFeatures))) {
    
    featureValues <- cellFeatures[[featureIndex]]
    cellFeatMST <- cellMST[[featureIndex]]
    
    cellMST[[featureIndex]] <- updateFeatureMarkov(featureValues, cellFeatMST)    
    
    
  }
  
  return(cellMST)
  
  
}


#cellFeatMST <- list(list(history, countHistory, nextData = list(list(nextValue, probability, countNext, totCount))))

featurePrevisionMarkov <- function(cellFeatureHistory, cellFeatMST, step_numb) {
  
  vj <- cellFeatureHistory[1]
  
  if(length(cellFeatureHistory) > 1) {
    vk <- cellFeatureHistory[2:length(cellFeatureHistory)]
  }
  
  #tot_jumps <- NULL
  curr_jumps <- list(list(curr_val = vj,
                          curr_hist = vk,
                          prob = 1,
                          step = 0))
  
  #search at i distance
  for(i in c(1:step_numb)) {
    
    next_jumps <- NULL
    
    for(jump in curr_jumps) {
      
      #history of next state
      nex_hist <- c(jump$curr_val, jump$curr_hist[-length(jump$curr_hist)])
      
      for(featHistory in cellFeatMST) {
        
        if(sameHistory(nex_hist, featHistory$history, exact = T)) {
          
          #predict next steps with probability
          for(value in featHistory$nextData) {
            
            next_jumps <- c(next_jumps, list(list(curr_val = value$nextValue,
                                                  curr_hist = nex_hist,
                                                  prob = jump$prob*value$probability,
                                                  step = i)))
            
          }
        }
      }
    }
    
    curr_jumps <- next_jumps
    #tot_jumps <- c(tot_jumps, next_jumps)
  }
  
  #Here aggregate data
  
  return(curr_jumps)
  
}

miniTestUpdate <- function() {
  
  CMST <- updateFeatureMarkov(c(4,3,2,1), NULL)
  CMST <- updateFeatureMarkov(c(5,4,3,2), CMST)
  CMST <- updateFeatureMarkov(c(6,5,4,3), CMST)
  CMST <- updateFeatureMarkov(c(6,7,9,2), CMST)
  CMST <- updateFeatureMarkov(c(4,3,2,1), CMST)
  CMST <- updateFeatureMarkov(c(5,3,2,1), CMST)
  CMST <- updateFeatureMarkov(c(5,3,2,1), CMST)
  CMST <- updateFeatureMarkov(c(7,4,3,2), CMST)
  CMST <- updateFeatureMarkov(c(5,4,3,2), CMST)
  CMST <- updateFeatureMarkov(c(5,4,3,2), CMST)
  
  return(CMST)
  
}

#returns updated cellFeatMST
updateFeatureMarkov <- function(cellFeatureHistory, cellFeatMST) {
  
  
  m <- 1
  b <- 4
  
  
  vj <- cellFeatureHistory[1]
  if(length(cellFeatureHistory) > 1) {
    vk <- cellFeatureHistory[2:length(cellFeatureHistory)]
  }
  
  totvj <- 1
  totnotvj <- 0
  totvk <- 0
  totnotvk <- 0
  
  #calculate total vj
  if(!is.null(cellFeatMST)) {
    for(featHistory in cellFeatMST) {
      if(!is.null(featHistory$nextData)) {
        for(value in featHistory$nextData) {
          if(sameValue(vj, value$nextValue)) {        
            totvj <- totvj + value$countNext
          }
        }
      }
    }
  }
  #check for matching record history to update correct line
  updatedHistory <- FALSE
  
  
  if(!is.null(cellFeatMST)) {
    for(findex in c(1:length(cellFeatMST))) {
      
      
      featHistory <- cellFeatMST[[findex]]
      
      
      if(sameHistory(vk, featHistory$history, exact = T)) {
        
        featHistory$countHistory <- featHistory$countHistory + 1
        totvk <- totvk + featHistory$countHistory
        
        #check for matching next value
        updatedNextValue <- FALSE
        
        if(!is.null(featHistory$nextData)) {
          for(vindex in c(1:length(featHistory$nextData))) {
            
            value <- featHistory$nextData[[vindex]]
            
            if(sameValue(vj, value$nextValue)) {
              
              value$countNext <- value$countNext + 1
              value$totCount <- value$totCount + 1
              updatedNextValue <- TRUE
              
            }
            else {
              totnotvj <- totnotvj + value$countNext
            }
            
            featHistory$nextData[[vindex]] <- value
          }
        }
        #add new next value
        if(updatedNextValue == FALSE) {
          
          newData <- list(list(nextValue = vj,
                               probability = 0,
                               countNext = 1,
                               totCount = totvj))
          
          featHistory$nextData <- c(featHistory$nextData, newData)
          
        }
        
        updatedHistory <- TRUE
        
      }
      else {
        #not vk entry
        totnotvk <- totnotvk + featHistory$countHistory
        
        #search for vj in not vk
        if(!is.null(featHistory$nextData)) {
          for(vindex in c(1:length(featHistory$nextData))) {
            value <- featHistory$nextData[[vindex]]
            if(sameValue(vj, value$nextValue)) {
              value$totCount <- value$totCount + 1
              
            }
            else {
              totnotvj <- totnotvj + value$countNext
              
            }
            featHistory$nextData[[vindex]] <- value
          }
        }
        
      }
      
      cellFeatMST[[findex]] <- featHistory
    }
  }
  
  #add new history
  if(updatedHistory == FALSE) {
    
    newData <- list(list(nextValue = vj,
                         probability = 0,
                         countNext = 1,
                         totCount = totvj))
    
    newHist <- list(list(history = vk,
                         countHistory = 1,
                         nextData = newData))
    
    totvk <- totvk + 1
    
    cellFeatMST <- c(cellFeatMST, newHist)
    
  }
  
  
  tothistories <- totvk + totnotvk
  totvalues <- totvj + totnotvj
  #updateProbabilities
  
  if(!is.null(cellFeatMST)) {
    for(findex in c(1:length(cellFeatMST))) {
      
      featHistory <- cellFeatMST[[findex]]
      
      Clog <- log(featHistory$countHistory / (tothistories))
      C <- (featHistory$countHistory / (tothistories))
      
      
      if(!is.null(featHistory$nextData)) {
        for(vindex in c(1:length(featHistory$nextData))) {
          
          value <- featHistory$nextData[[vindex]]
          
          Blog <- log(value$totCount / totvalues)
          B <- (value$totCount / totvalues)
          
          Alog <- log((value$countNext +(m/b))/ (value$totCount + m))      
          A <- ((value$countNext +(m/b))/ (value$totCount + m))   
          
          value$probability <- A*B/C
          #logprob
          #value$probability <- (Alog + Blog) - Clog
          
          featHistory$nextData[[vindex]] <- value
          
        }
      }
      
      cellFeatMST[[findex]] <- featHistory
      
    }
    
  }
  
  
  return(cellFeatMST)
}


sameValue <- function(val1, val2) {
  
  return(val1 == val2)
  
}

#exact: if T only returns perfect matches. If F NA is considered wild card
#Use exact = T when building markov chain, may use exact = F when searching it
sameHistory <- function(h1, h2, exact = TRUE) {
  
  res <- all(h1 == h2, na.rm = !exact)
  
  if(is.na(res)) {
    return(FALSE)
  }
  else {
    return(res)
  }
  
  
}


###################### Fast Subset Scan ###################################

testFSS <- function() {
  
  node_number = 5
  feat_number = 4
  feat_memory = 12
  
  learning_iterations = 4
  
  period = 24
  
  feat_combinations <- powerset(c(1:feat_number))
  
  periodHistory <- rep(list(rep(list(matrix(nrow = feat_number, ncol = feat_memory)), period)), node_number)
  residAvg <- rep(list(vector("numeric", length = feat_number)), node_number)
  residVar <- rep(list(vector("numeric", length = feat_number)), node_number)
  
  monitorResiduals <- rep(list(NULL), node_number)
  
  
  nodeList <- createTestNodeList(node_number, feat_number, period, period_number = 50)
  
  #add anomaly
  nodeList[[1]][2,1000] <- nodeList[[1]][2,1000] * 0.4
  nodeList[[3]][2,1000] <- nodeList[[3]][2,1000] * 0.4
  nodeList[[4]][2,1000] <- nodeList[[4]][2,1000] * 0.4
  
  for(time in c(1:(ncol(nodeList[[1]])))) {
    
    currentPeriod <- ((time-1) %% period) + 1
    print(paste("Time ",time,", hour ",currentPeriod))
    
    if(time == period * feat_memory + learning_iterations * period + 1) {
      print("Start node scoring")
    }
    
    
    
    currentValues <- rep(list(NULL), node_number)
    baseValues <- rep(list(NULL), node_number)
    
    for(feat_comb in feat_combinations) {
      
      node_scores <- NULL
      
      for(node in c(1:node_number)) {
        
        
        currentPeriodHistory <- matrix(periodHistory[[node]][[currentPeriod]][feat_comb,], nrow = length(feat_comb))
        currentFeatures <- nodeList[[node]][feat_comb,time]
        currentResiduals <- NULL
        
        
        #Begin tests only after learning period
        
        if(time > period * feat_memory + learning_iterations * period) {
          
          
          muVals <- rowMeans(currentPeriodHistory)
          sigVals <- rowSums((currentPeriodHistory - muVals)^2)/(dim(currentPeriodHistory)[2] - 1)
          ranl <- apply(currentPeriodHistory, 1, range)
          ranVals <- ranl[2,] - ranl[1,]
          
          currentResiduals <- (currentFeatures - muVals) / ranVals
          
          
          
          node_multiv_score <- node_scoring(currentPeriodHistory, currentFeatures, residAvg[[node]][feat_comb], residVar[[node]][feat_comb])
          
          node_scores <- c(node_scores, node_multiv_score)
        }
        else if(time > period * feat_memory) {
          
          #learn values for residual modeling
          
          muVals <- rowMeans(currentPeriodHistory)
          sigVals <- rowSums((currentPeriodHistory - muVals)^2)/(dim(currentPeriodHistory)[2] - 1)
          ranl <- apply(currentPeriodHistory, 1, range)
          ranVals <- ranl[2,] - ranl[1,]
          
          currentResiduals <- (currentFeatures - muVals) / ranVals
          
          
          time_iteration <- (time - period*feat_memory)
          
          lastAvg <- residAvg[[node]][feat_comb]
          
          residAvg[[node]][feat_comb] <- (residAvg[[node]][feat_comb] * (time_iteration - 1) + 
                                            currentResiduals) / time_iteration
          
          residVar[[node]][feat_comb] <- (((time_iteration - 1) * residVar[[node]][feat_comb])
                                          + (currentResiduals - lastAvg) * (currentResiduals - residAvg[[node]][feat_comb])) / time_iteration        
          
        }
        
        if(!is.null(currentResiduals)) {
          
          #monitorResiduals[[node]] <- cbind(monitorResiduals[[node]], currentResiduals)
          currentValues[[node]] <- currentResiduals
          
        }
        
        baseValues[[node]] <- cbind(residAvg[[node]][feat_comb], residVar[[node]][feat_comb])
        
      }
      
      
      #cluster of a node with these neighbors
      
      if(time > period * feat_memory + period * learning_iterations) {
        
        
        cluster_scoring(priorities = node_scores, 
                        current = currentValues, 
                        base = baseValues)
        
        
      }
      
      
      
    }
    
    
    for ( node in c(1:node_number)) {
      
      currentPeriodHistory <- shiftRightMatrix(currentPeriodHistory, 1)
      currentPeriodHistory[, 1] <- nodeList[[node]][,time]
      periodHistory[[node]][[currentPeriod]] <- currentPeriodHistory
      
    }
    
    
  }
  
  #plot(monitorResiduals[[3]][2,], type = "o")
  
}



createTestNodeList <- function(node_number, feat_number, period = 24, period_number = 8) {
  
  nodeList <- NULL
  #5 nodes, each with 4 features, for 192 time slot
  for(i in c(1:node_number)) {
    
    node <- createTestNode(mean = 0, range = 2, noise_level = 0.1, 
                           feature_number = feat_number, period, period_number)
    
    nodeList <- c(nodeList, list(node))
    
  }
  return (nodeList)
  
}

createTestNode <- function(mean, range, noise_level, feature_number, period, period_number) {
  
  node <- NULL
  
  for(i in c(1:feature_number)) {
    
    sinusoid <- sin(pi*(0:((period*period_number) - 1))/(period/2)) * (range / 2) + mean +
      #  runif(period*period_number, -noise_level, noise_level)
      rnorm(period*period_number, mean = 0, sd = noise_level)
    
    node <- rbind(node, sinusoid)
    
  }
  
  return(node)
  
  
}

#FGSS functions 

#for multivariate: pass the correct subset of features
#aggregated by feature
#period history: for each j feature, remember m past elements from previous periods
node_scoring <- function(currentPeriodHistory, currentFeatures, residAvg, residVar) {
  
  
  mu <- rowMeans(currentPeriodHistory)  
  sig2 <- rowSums((currentPeriodHistory - mu)^2)/(dim(currentPeriodHistory)[2] - 1)
  
  ranl <- apply(currentPeriodHistory, 1, range)
  ran <- ranl[2,] - ranl[1,]
  
  
  
  residuals <- (currentFeatures - mu) / ran
  
  Cim <- (residuals - residAvg) ^ 2 / residVar
  Bim <- residVar / residVar
  
  
  #aggregate Cim and Bim
  Ci <- sum(Cim)
  Bi <- sum(Bim)
  
  
  priority <- Ci
  
  return(priority)
  
  
}


#priorities: priority (node score) for each n node.
#current: list of n, values of current features for node.
#base: list of n, matrix of j rows of column mu and sig2 

cluster_scoring <- function(priorities, current, base, anomThreshold = 5) {
  
  #order priorities
  prior_order <- sort(priorities, decreasing = TRUE,index.return = T)$ix
  
  subset <- NULL
  C <- rep( 0  , length(prior_order))
  B <- rep( 0  , length(prior_order))
  
  P <- rep( 0  , length(prior_order))
  
  
  #in order, calculate F(S)
  for (size in c(1:length(prior_order))) {
    
    node <- prior_order[size]
    
    residuals <- current[[node]]
    resAvgs <- base[[node]][,1]
    resVars <- base[[node]][,2]
    
    Cim <- (residuals - resAvgs) ^ 2 / resVars
    Bim <- resVars / resVars
    
    Ci <- sum(Cim)
    Bi <- sum(Bim)
    
    C[size:length(C)] <- C[size:length(C)] + Ci
    B[size:length(B)] <- B[size:length(B)] + Bi
    
    P[size:length(P)] <- P[size:length(P)] + priorities[node]
  }
  
  #final_scores <- (C - B)^2 / (2*B)
  
  final_scores <- 0.5 * ( B * log( B / C ) + C - B )
  
  string <- NULL
  
  #   #Full results
  #   for(size in c(1:length(prior_order))) {
  #     
  #     string <- paste(string," ", prior_order[size])
  #     
  #     print(paste("Set ",size,", nodes (",string,"), score ",final_scores[size]))
  #     
  #     
  #   }
  #    
  
  
  #Best result
  
  bestSize <- 0
  bestScore <- 0
  
  for(size in c(1:length(prior_order))) {
    if(final_scores[size] > bestScore) {
      bestScore <- final_scores[size]
      bestSize <- size  
    }
  }
  for(size in c(1:bestSize)) {
    string <- paste(string," ", prior_order[size])
  }
  
  if(bestScore > anomThreshold) {
    print(paste("Set ",size,", nodes (",string,"), score ",bestScore))
  }
  
  
}

powerset <- function(s){
  
  len <- length(s)
  l <- vector(mode="list",length=2^len) ; 
  l[[1]]=numeric()
  counter <- 1
  for(x in 1:length(s)) {
    for(subset in 1:counter){
      counter <- counter+1
      l[[counter]] <- c(l[[subset]],s[x])
    }
  }
  l[[1]] <- NULL
  return(l)
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
  #   
  #     sums <- (nrow(myData)-1)*sum(apply(myData,1,var))
  #   AIC <- vector("numeric", length = nrow(myData))
  #   fitty <- NULL
  #   
  #   for ( i in (1:(nrow(distMatrix) - 1))) {
  #     
  #    myFit <- kmeans(distMatrix, centers = i)
  #     
  #     AIC[i] <- FOM(myFit)
  #     fitty <- c(fitty, list(myFit))
  #     
  #     sums[i] <- sum(kmeans(as.matrix(distMatrix), centers = i)$withinss)
  #     
  #   }  
  # 
  # 
  #   #auto search elbow point
  # 
  #   
  #   v = -diff(AIC)
  #   nv = length(v)
  #   fom = v[1:(nv-1)]/v[2:nv]
  #   nclus = which.max(fom)+1
  #   cat("The apparent number of clusters is: ",nclus,"\n")
  # #   print("Cluster membership: ")
  # #   print(paste(fitty[[nclus]]$cluster))
  #   
  #   plot(AIC, type = "o")
  #   points(nclus,AIC[nclus],col=2,pch=20,cex=2)
  #   
  
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


###################### Weighting Fun ############


pgam <- function(shape, scaleFactor) {
  
  i <- 1:100
  i <- i / 100
  
  plot(i,i, type = "l")
  points(i, pgamma(i, shape = shape, scale = 1 / (shape * scaleFactor)), col = "red", type = "l")
  
}

qgam <- function(shape, scaleFactor, scale = NULL) {
  
  i <- 1:100
  i <- i / 100
  
  plot(i,i, type = "l")
  if(is.null(scale)) {
    points(i, qgamma(i, shape = shape, scale = 1 / (shape * scaleFactor)), col = "red", type = "l")
    
  }
  else {
    points(i, qgamma(i, shape = shape, scale = scale), col = "red", type = "l")
    
  }
}

pwr <- function(x, exp) {
  
  return (sign(x) * abs(x) ^ exp)
  
  #return (sign(x) * exp ^ abs(x))
  
  
}

ppower <- function(confid, weight) {
  
  i <- 1:10000
  i <- i / 10000
  
  
  correctConfid <- confid
  confidIndex <- which(i == confid)
  
  myPower <- pwr( (i - correctConfid), 1 / weight) + confid
  #myPower2 <- pwr(i[10000:1] - correctConfid, 1 / weight) + confid
  dist <- myPower - i
  #dist2 <- myPower2 - i
  #   #myPower  <- (myPower * 1/weight)  + confid
  
  #myPower <- myPower - dist * (abs(i - confid)/confid)
  myT <- myPower - (dist * (abs(i-confid)/confid)^(1/weight))
  myT2 <- myPower - (dist * (abs(i-confid)/(1 - confid))^(1/weight))
  #myT3 <- myPower2 - (dist2 * (abs(i[10000:1]-confid)/confid))
  #myT4 <- myPower2 - (dist2 * (abs(i[10000:1]-confid)/(1 - confid)))
  #   myPower <- ( myPower - min(myPower))/ (max(myPower) - min(myPower))
  
  #   
  #   confidIndex <- which(myT == confid)
  myPower[1:(confidIndex - 1)] <- myT[1:(confidIndex - 1)]
  myPower[(confidIndex + 1):length(myPower)] <- myT2[(confidIndex + 1):length(myPower)]
  
  #   myPower2[1:(confidIndex - 1)] <- myT3[1:(confidIndex - 1)]
  #   myPower2[(confidIndex + 1):length(myPower2)] <- myT4[(confidIndex + 1):length(myPower)]
  # 
  #   
  #   #myPower[1:confidIndex] <- myPower[1:confidIndex] - ((confidIndex - c(1:confidIndex)) / confidIndex)  * myPower[1]
  #myPower[confidIndex:length(myPower)] <- myPower[confidIndex:length(myPower)] + (abs(confidIndex - c(confidIndex:length(myPower)))/ abs(length(myPower) - confidIndex)) * (1 - myPower[length(myPower)])
  
  
  plot(i,i, type = "l")
  #points(i, (pwr((i - confid) , 1 / weight) + confid), type = "l", col = "red")
  #points(i, myPower, type = "l", col = "red")
  
  points(i, myPower, col = "green", type = "l")
  #points(i, myPower2, col = "purple", type = "l")
  #points(i, myT, col = "blue", type = "l")
  #points(i, myT2, col = "purple",type = "l")
  
  
  
  
}

sigmoid <- function(confid, k) {
  
  i <- 1:10000
  i <- i / 10000
  
  pin <- ( 1 / (1 + exp((k * -(i - confid)))  ))
  pin <- (pin - min(pin)) / (max(pin) - min(pin))
  
  reflex <- i + (i - pin)
  
  
  plot(i,i,type = "l")
  points(i, pin, col = "red", type = "l")
  
  
  
}

logit <- function(confid, k) {
  
  i <- 1:10000
  i <- i / 10000
  
  pin <- - log( (1 / k*(i - confid) - 1))
  pin <- (pin - min(pin)) / (max(pin) - min(pin))
  
  reflex <- i + (i - pin)
  
  
  plot(i,i,type = "l")
  points(i, pin, col = "red", type = "l")
  
  
  
}


###################### PSF + Motifs #############


#Single time series (or multi? in that case on which subset of neighbors, and how to incorporate weights?)
#Input: a single time series (do not deseason), a vector
psf_label_set <- function(training_set, min_k = 1, max_k = 20) {
  
  lab_set <- NULL
  lev_set <- NULL
  
  for(row in (1:nrow(training_set))) {
    lab_serie <- psf_label_series(training_set[row,], min_k, max_k)
    lev <- levels(lab_serie)
    
    lab_set <- rbind(lab_set, lab_serie)
    lev_set <- c(lev_set, list(lev))
    
  }
  
  return(list(labels = lab_set,
              levels = lev_set))
}

psf_label_series <- function(training_series, min_k = 1, max_k = 20) {
  
  library(fpc) # for pamk
  
  #pamk to find best k labels to describe (original: day wrt to year; mine: hour\ set of hours)
  
  if(min_k < max_k) {
    pamcl <- pamk(training_series, krange = c(min_k:max_k))
    labeled_training_series <- pamcl[[1]]$clustering
    bestK <- pamcl$nc
  }
  
  else {
    bestK <- min_k
  }
  
  
  
  factorial <- cut(training_series, bestK)
  labeled_training_series <- factorial
  
  return(labeled_training_series)
  
}



#Single time series \ multi
#Input: A set of NEIGHBORING (this is to be done BEFORE psf_forecast) training_series, 
# their labeled counterpart (via psf_label_set), the window_size to be used, 
# the index of the day to be forecast, weights according to similarity of neighbors. 
# The first series in training_set is considered the calculating node.
psf_forecast<- function(training_set, label_set, window_size, current_day, offset = 1, weight_set = NULL, window_change = T) {
  
  main_dist_threshold <- 1
  neigh_dist_threshold <- 2
  
  #PSF
  
  #assumes current_day >= window_size + offset
  
  #find in training_set, windows with same labels as window before current_day
  ##check if windows of nearest neighbors are also similar
  
  confirmed_match <- NULL
  confirmed_match_number <- 0
  weight_match <- NULL
  
  perfect_match <- 1
  close_match <- 0.8
  far_match <- 0.5
  
  while(confirmed_match_number <= 0 && window_size >= 1) {
    
    windowed_set <- as.matrix(label_set[, c((current_day - window_size - offset + 1) : (current_day - offset))])
    
    
    for(t in c((window_size + offset):ncol(label_set))) {
      if(t != current_day) {
        
        #How many categorical labels are different?
        main_distance <- length(which(windowed_set[1,] != label_set[1,((t - window_size - offset + 1) : (t - offset))]))
        neigh_distance <- length(which(windowed_set[-1,] != label_set[-1,((t - window_size - offset + 1) : (t - offset))]))
        
        if(main_distance <= main_dist_threshold) {
          
          #found matching self
          
          confirmed_match <- c(confirmed_match, training_set[1,t])
          confirmed_match_number <- confirmed_match_number + 1 
          
          if(neigh_distance <= neigh_dist_threshold) {
            
            #matches both own history and neigh history
            weight_match <- c(weight_match, perfect_match)
            
          }
          else {
            
            #matches own history, not neighborhood
            weight_match <- c(weight_match, close_match)
            
          }
          
        }
        else if(neigh_distance <= neigh_dist_threshold) {
          
          confirmed_match <- c(confirmed_match, training_set[1,t])
          confirmed_match_number <- confirmed_match_number + 1 
          
          #matches neigh history, not own... outlier?
          weight_match <- c(weight_match, far_match)
          
        }
      }
    }
    
    if(confirmed_match_number <= 0 && window_change == T) {
      #retry with smaller window
      window_size <- window_size - 1
      
    }
    else if(window_change == F) {
      break
    }
    
    
  }
  
  #forecast average between all found values
  ##weight average according to similarity of context too (nearest neighbor window)
  
  fc <- weighted.mean(confirmed_match, weight_match)
  
  return(fc)
  
}


label_with_interval <- function(training_set, interval_set) {
  
  
  lab_set <- NULL
  
  for(i in (1:nrow(training_set))) {
    
    intervals <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", interval_set[[i]]) ),
                       upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", interval_set[[i]]) ))
    
    new_set <- vector("numeric", length = length(training_set[i,]))
    new_set <- rep( nrow(intervals) + 1 , length(training_set[i,]))
    
    for(j in (nrow(intervals) : 1)) {
      #find indices
      indices <- which(training_set[i,] <= intervals[j,2])
      #set the label if not already present
      new_set[indices] <- min(new_set[indices], j)
      
    }
    
    lab_set <- rbind(lab_set, new_set)
    
  }
  
  return(lab_set)
  
  
}

bestWindowDiscovery <- function(training_set, level_set, n, offset = 1, min_window = 1, max_window = 24, minK = 1, maxK = 24, return_error = FALSE) {
  
  #training_set subdivsion, for cross_validation
  
  #Because of the time relationship, must split the set in meaningful chunks
  cross_validation_sets <- NULL
  #note that the smaller the resulting width after dividing the set, the smaller the
  #window size can be!    max_window_size <= width - offset   !!! 
  
  width <- floor(ncol(training_set) / n)
  
  for(i in (1:n)) {
    
    cross_validation_sets <- c(cross_validation_sets,
                               list(training_set[,(((i-1)*width + 1) : (i*width))]))
    
  }
  #OR signal error
  max_window_size <- max(min(max_window, width - offset), min_window)
  maxK <- min(maxK, width - 1)
  
  best_window_error <- Inf
  best_window_size <- 0
  
  #Find best window size
  for(window_size in c(min_window:max_window_size)) {
    
    expected_fold_errors <- NULL
    
    for(n_fold in (1:n)) {
      
      n_set <- cross_validation_sets[[n_fold]]
      if(is.null(level_set)){
        lab_set <- psf_label_set(n_set,minK,maxK)$labels
        
      }
      else{
        lab_set <- label_with_interval(n_set, level_set)
        
      }
      time_errors <- NULL
      for(t in ((window_size + offset):ncol(n_set))) {
        
        #forecast and calculate error
        forec <- psf_forecast(n_set, lab_set, window_size, current_day = t, offset = offset, window_change = F)
        #calc error of day in remaining month
        
        time_errors <- c(time_errors, abs(forec - n_set[1,t]))
        
      }
      #calculate error with this window size for this fold
      no_forecast_fraction <- length(which(is.nan(time_errors))) / length(time_errors)
      average_error <- mean(time_errors, na.rm = T)
      
      
      if(is.nan(average_error)) {
        experr <- no_forecast_fraction
        
      }
      else{
        experr <- sqrt((average_error)^2 + (no_forecast_fraction)^2)
      }
      
      expected_fold_errors <- c(expected_fold_errors, experr)
      
    }
    #calculate expected error for this window size
    expected_window_error <- mean(expected_fold_errors)
    
    #find best window size, minimizing expected_window_error
    if(expected_window_error < best_window_error) {
      best_window_size <- window_size
      best_window_error <- expected_window_error
    }
    
    
  }
  
  if(return_error == FALSE) {
    return(best_window_size)
  }
  else {
    
    return(list(win_size = best_window_size,
                exp_error = best_window_error))
    
  }
}




motif_discovery <- function(training_set, n, offset = 1, max_window = 24, min_k = 1, max_k = 24) {
  
  
  
  motif_set <- NULL
  candidate_set <- NULL
  outlier_set <- NULL
  
  print("Finding label set...")
  lab_set <- psf_label_set(training_set, min_k, max_k)$labels
  print("Finding best window...")
  window_disc <- bestWindowDiscovery(training_set, n, offset, max_window, min_k, max_k, return_error = TRUE)
  
  window_size <- window_disc$win_size
  exp_error <- window_disc$exp_error
  
  print(paste("Best window ",window_size))
  
  time_errors <- rep(0, window_size + offset - 1)
  print("Finding candidate set...")
  #for each time slot, calculate error and find candidate_set
  for(time in ((window_size + offset): ncol(training_set))) {
    
    forec <- psf_forecast(training_set, lab_set, window_size, current_day = time, offset = offset)
    error <- abs(forec - training_set[1,time])
    
    time_errors <- c(time_errors, error)
    
    #if error more than average, add time to candidate set
    if(error > exp_error) {
      candidate_set <- c(candidate_set, time)
    }
    
    
  }
  
  #deal with NaN ..?
  #NaN = not predicted = max error, unique motifs
  nan_candidate <- which(is.nan(time_errors))
  
  #cluster candidate_set over error, use highest priority cluster the ones with higher deviation from forecast (error)
  print("Clustering candidates...")
  cluster_candidate <- kmeans(time_errors[candidate_set], centers = 3)
  sorted_cluster_means <- sort(cluster_candidate$centers, decreasing = TRUE, index.return = TRUE)$ix
  candidate_set_h <- candidate_set[which(cluster_candidate$cluster == sorted_cluster_means[1])]
  candidate_set_m <- candidate_set[which(cluster_candidate$cluster == sorted_cluster_means[2])]
  candidate_set_l <- candidate_set[which(cluster_candidate$cluster == sorted_cluster_means[3])]
  
  
  main_dist_threshold <- 0
  neigh_dist_threshold <- 0
  
  
  
  print("Isolating meaningful motifs...")
  for(cand in (candidate_set_h)) {
    
    #find motifs (windows) before candidate_set.
    ##find windows of nn too
    candidate_motif <- as.matrix(lab_set[, (cand - offset - window_size +1):(cand - offset)])
    
    
    matches <- 1
    confirmed <- TRUE
    
    for(t in c((window_size + offset):ncol(lab_set))) {
      
      #How many categorical labels are different?
      main_distance <- length(which(candidate_motif[1,] != lab_set[1,((t - window_size - offset + 1) : (t - offset))]))
      neigh_distance <- length(which(candidate_motif[-1,] != lab_set[-1,((t - window_size - offset + 1) : (t - offset))]))
      
      if(main_distance <= main_dist_threshold) {
        
        #found match
        #If the motif only appear before candidates, consider it a valid motif (with given support), if not discard it
        ##check motif of nn too
        
        
        if(neigh_distance <= neigh_dist_threshold) {
          
          if(!(t %in% candidate_set)) {
            #motif appears before normal data
            confirmed <- FALSE
            
          }
          else {
            #perfect match on another candidate
            #increase support
            matches <- matches + 1
          }
          
        }
        else {
          
          #matches own history, not neighborhood        
          
          if(!(t %in% candidate_set)) {
            #motif appears before normal data
            confirmed <- FALSE
            
          }
          else {
            #close match on another candidate
            #increase support
            matches <- matches + 1
          }
          
        }
        
      }
      else if(neigh_distance <= neigh_dist_threshold) {
        
        
        #matches neigh history, not own... outlier?            
        if(!(t %in% candidate_set)) {
          #motif appears before normal data
          confirmed <- FALSE
          
        }
        else {
          #far match on another candidate
          #increase support
          matches <- matches + 1
        }
        
      }
    }
    
    if(confirmed == TRUE) {
      
      outlier_set <- c(outlier_set, cand)
      motif_set <- c(motif_set, list(candidate_motif))
      
    }
  }
  
  #   
  #   plot(lab_set)
  #   points(outlier_set, lab_set[outlier_set], col = "red")
  #   points()
  #   
  
  return(motif_set)
  
  
}




############## Sinusoids ##############################

#build data sets ( sinusoids ) of various kinds for testing of noisy and single anomaly detection

createTestData1 <- function(period = 24, period_number = 7) {
  
  #normal sinusoid with chosen period, for a certain period_number, with lowest value 0.5, highest 2.5 
  sinusoid1 <- sin(pi*(0:((period*period_number) - 1))/(period/2)) + 1.5
  
  #robustness to noise
  
  #sinusoid with low noise
  sinusoid2 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  #sinusoid with medium noise
  sinusoid3 <- sinusoid1 + runif(period*period_number, -0.4, 0.4)
  #sinusoid with high noise
  sinusoid4 <- sinusoid1 + runif(period*period_number, -0.5, 0.8)
  
  
  #single anomaly detection
  
  #single point, low noise, high peak of data, small anomaly score
  sinusoid5 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid5[32] <- sinusoid5[32] * c(0.6)
  #single point, low noise, high peak of data, high anomaly score
  sinusoid6 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid6[54] <- sinusoid6[54] * c(0.3)
  
  #single point, low noise, rising slope of data, small anomaly score
  sinusoid11 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid11[51] <- sinusoid11[51] * c(0.6)
  #single point, low noise, rising slope of data, high anomaly score
  sinusoid12 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid12[75] <- sinusoid12[75] * c(0.3)
  
  #single point, low noise, downward slope of data, small anomaly score
  sinusoid13 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid13[83] <- sinusoid13[83] * c(0.6)
  #single point, low noise, downward slope of data, high anomaly score
  sinusoid14 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid14[107] <- sinusoid14[107] * c(0.3)
  
  #single point, low noise, low valley of data, small anomaly score
  sinusoid7 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid7[87] <- sinusoid7[87] * c(0.6)
  #single point, low noise, low valley of data, high anomaly score
  sinusoid8 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid8[113] <- sinusoid8[113] * c(0.3)
  
  
  #single point, high noise, high peak of data, high anomaly score
  sinusoid9 <- sinusoid1 + runif(period*period_number, -0.5, 0.5)
  sinusoid9[103] <- sinusoid9[103] * c(0.5)
  #single point, high noise, low valley of data, high anomaly score
  sinusoid10 <- sinusoid1 + runif(period*period_number, -0.5, 0.5)
  sinusoid10[135] <- sinusoid10[135] * c(0.5)
  
  
  clust <- rbind(sinusoid1,
                 sinusoid2,
                 sinusoid3,
                 sinusoid4,
                 sinusoid5,
                 sinusoid6,
                 sinusoid7,
                 sinusoid8,
                 sinusoid9,
                 sinusoid10,
                 sinusoid11,
                 sinusoid12,
                 sinusoid13,
                 sinusoid14)
  
  return(clust)
  
  
}


createTestData2 <- function(period = 24, period_number = 7) {
  
  
  #normal sinusoid with chosen period, for a certain period_number, with lowest value 0.5, highest 2.5 
  sinusoid1 <- sin(pi*(0:((period*period_number) - 1))/(period/2)) + 1.5
  
  
  #anomaly length 5, fast down, fast up
  sinusoid15 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid15[83:87] <- sinusoid15[83:87]*c(0.4, 0.2, 0.2, 0.2, 0.8)
  #anomaly length 5, fast down, slow up
  sinusoid16 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid16[83:87] <- sinusoid16[83:87]*c(0.4, 0.55, 0.65, 0.7, 0.85)
  #anomaly length 5, slow down, fast up
  sinusoid17 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid17[83:87] <- sinusoid17[83:87]*c(0.8, 0.65, 0.5, 0.3, 0.8)
  #anomaly length 5, slow down, slow up
  sinusoid18 <- sinusoid1 + runif(period*period_number, -0.3, 0.3)
  sinusoid18[83:87] <- sinusoid18[83:87]*c(0.75, 0.6, 0.4, 0.55, 0.7)
  
  clust <- rbind(sinusoid1,
                 sinusoid15,
                 sinusoid16,
                 sinusoid17,
                 sinusoid18)
  
  return(clust)
  
  
  
}


#calculates probability that a single value is a good fit to the distribution assumed 
#from the rest of the nodes

fitnessCheck <- function(historic, newValue, confidence, plot = FALSE) {
  
  library(forecast)
  
  #fit an arima to the historical data
  aa <- auto.arima(historic)
  ff <- forecast(aa, h = 1, level = confidence)
  
  if(plot == TRUE) {
    plot(ff)
    lines(rep(newValue, length(historic) + 1), col = "red")
  }
  if(newValue >= ff$lower[1] && newValue <= ff$upper[1]) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
  
}


############## TEST SET GENERATION ###########

#Generates a single time series with given parameters
generateTimeSeries <- function(mean = 0, range = 1, periodicity = 24, noise = range/5, end = periodicity * 30, start_shift = 0, type = "sinusoid") {
  
  series <- NULL
  
  if(type == "sinusoid") {
    series <- sin(pi*(c((0:(end-1)) + start_shift)/(periodicity/2)) ) * range + mean + runif(end, -noise, +noise)
  }
  return(series)
  
  
}

#Generates a matrix of n time series related to same feature
generateFeatureSets <- function(location_number, mean = 0, range = 1, periodicity = 24, noise = range/5, end = periodicity * 30, start_shift = 0, type = "sinusoid") {
  
  set <- NULL
  
  for(i in (1:location_number)) {
    
    set <- rbind(set, generateTimeSeries(mean, range, periodicity, noise, end, start_shift, type))
    
  }
  
  return(set)
  
  
}


#adds a random noise to a series, between start and end indices
#dir values: "two" adds noise in both directions; "up" and "down" adds only positive or negative noise
series.addNoise <- function(series, noise = ( max(series) - min(series) ) / 10, 
                            start = 1, end = length(series), dir = "two") {
  
  min_noise <- -noise
  max_noise <- noise
  
  if(dir == "up") {
    min_noise  <- 0
  }
  if(dir == "down") {
    max_noise <- 0    
  }
  
  series[start:end] <- series[start:end] + runif(length(start:end), min_noise, max_noise)
  
  
  return(series)
  
}

#generates random coordinates
generateCoordinates <- function(loc_number, min_x = 0, max_x = 1, min_y = 0, max_y = 1) {
  
  coord <- NULL
  
  for(i in (1:loc_number)) {
    
    newcoord <- c(runif(1, min_x, max_x), runif(1, min_y, max_y))
    
    coord <- rbind(coord, newcoord)
    
    
  }
  
  return(coord)
  
}

#adds a random noise to all series, between start and end indices, within
# a range from a center

set.addNoiseArea <- function(set, coordinates, center, range, noise = NULL, start = 1, end = ncol(set), dir = "two", proportional = F) {
  
  #find node_numbers of nodes within range to center
  distances <- dist(rbind(center,coordinates))[1:nrow(coordinates)]
  node_numbers <- which(distances < range)
  
  #add noise to selected nodes
  for(n in node_numbers) {
    if(is.null(noise)) {
      noise  <- (max(set[n,]) - min(set[n,])) / 10
    }
    if(proportional == T) {
      noise <- noise * (1 - (distances[n] / range))
      
    }
    set[n,] <- series.addNoise(set[n,], noise, start, end, dir)
  }
  
  print(paste("Added noise to nodes: ", node_numbers))
  return(set)
  
}


