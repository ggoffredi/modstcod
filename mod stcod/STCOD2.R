
#cluster = a table with rows of nodes in a cluster
STCOD_Cluster <- function(cluster, period = 24, watch_node = 5, neigh_number = 5) {
  
  state_mem <- 5
  
  genAnomaly <- FALSE
  
  anomalyPredictionOn <- rep(FALSE, nrow(cluster))
  predictMode <- rep(FALSE,nrow(cluster))
  predictedAnomaly <- rep(FALSE,nrow(cluster))
  predictNext <- rep(list(),nrow(cluster))
  
  cell_threshold <- 2.5
  diff_threshold <- 2
  vote_threshold <- 1
  similarity_threshold <- 0.5
  
  cell_update_rate_OK <- 0.6
  cell_update_rate_OUT <- 0.25
  
  cluster_update_rate <- 0.3
  
  #for simplicity test without na
  cluster[is.na(cluster)] <- 0
  
  #following data is stored in list instead of matrix to simulate cell distributed behavior
  #element in a list[[1]] represents knowledge of node 1
  #this is redundant in this algorithm, as cluster_season_sample and prev_state should be the same
  
  #seasonality data
  cluster_seasonality <- rep(list(vector("numeric", length = period)), nrow(cluster))
  cell_seasonality <- rep(list(vector("numeric", length = period)), nrow(cluster))
  
  cell_season_diff <- rep(list(vector("numeric", length = period)), nrow(cluster))
  cell_residual <- rep(list(vector("numeric", length = period)), nrow(cluster))
  cell_val_std <- rep(list(vector("numeric", length = period)), nrow(cluster))
  cell_acalc <- rep(list(vector("numeric", length = period)), nrow(cluster))
  cell_season_lon_diff <- rep(list(vector("numeric", length = period)), nrow(cluster))
  
  #state memory of cell and neighbors
  #prev_state[[node]][1] indicates current state (after phase 1), [2] indicates state in previous time etc...
  prev_state <- rep(list(matrix(, nrow = nrow(cluster), ncol = state_mem)), nrow(cluster))
  prev_neigh_value <- rep(list(matrix(, nrow = nrow(cluster), ncol = state_mem)), nrow(cluster))
  neigh_weight <- rep(list(vector("numeric", length = nrow(cluster))), nrow(cluster))
  neigh_index <- rep(list(vector("numeric", length = nrow(cluster))), nrow(cluster))

  
  anomaly_structure <- rep(list(NULL), nrow(cluster))
  anomaly_unusedId <- 1
  
  
  
  #initialize weights and neighbors
  for(node in  (1:nrow(cluster))) {
    
    
    neigh_weight[[node]][] <- 1
    
    library(FNN)
    #should use spatial coord
    spnn <- get.knn(cluster, neigh_number, "kd_tree")
    nn.index <- spnn[[1]]
    
    neigh_index[[node]] <- nn.index[node,]
#     
#     
#     neigh_weight[[node]][] <- 0
#     neigh_weight[[node]][neigh_index[[node]]] <- 1
    
  }
  
  
  for(time in (1:ncol(cluster))) {
    
    
    
    #phase 1 - read and self detection
    
    #START NODE OPERATION
    for(node in (1:nrow(cluster))) {
      
      
      #check time slot
      curr_period <- ((time-1) %% period) + 1
      
      #read current and previous value (should be stored)
      curr_value <- cluster[node, time]
      if(time > 1) {
        prev_value <- cluster[node,time - 1]
      }
      else {
        prev_value <- 0
      }
      
      #get cell seasonality
      curr_value_season <- cell_seasonality[[node]][curr_period]
      if(curr_period == 1) {
        prev_value_season <- cell_seasonality[[node]][period]
        prev_season_diff <- cell_season_diff[[node]][period]
        prev_season_lon_diff <- cell_season_lon_diff[[node]][period]
        prev_residual <- cell_residual[[node]][period]
      }
      else {
        prev_value_season <- cell_seasonality[[node]][curr_period - 1]
        prev_season_diff <- cell_season_diff[[node]][curr_period - 1]
        prev_season_lon_diff <- cell_season_lon_diff[[node]][curr_period - 1]
        prev_residual <- cell_residual[[node]][curr_period - 1]
      }
      
      season_avg <- mean(cell_seasonality[[node]], na.rm = T)
      season_rng <- (range(cell_seasonality[[node]])[2] - 
                       range(cell_seasonality[[node]])[1])
      
      #if no season value exists, skip any check, consider node OK (learning phase)
      if(time < period + 1) {
        
        
        #update current state, shift previous for this node in its list
        prev_state[[node]][node,] <- shiftRight(prev_state[[node]][node,] , 1)
        prev_state[[node]][node,1] <- "OK"
        
        #initialize cell seasonality as same as the first period of values
        cell_seasonality[[node]][curr_period] <- cluster[node,time]
      }
      #else check seasonality (wrt cell past seasons)
      else {
        
        if(is.na(prev_season_diff)){
          prev_season_diff <- 0
        }
        
        #standardize comparison values from -1 to +1
        
        if(season_rng != 0) {
          curr_value_std <- ((curr_value - season_avg)/(season_rng/2))
          prev_value_std <- ((prev_value - season_avg)/(season_rng/2))
          curr_value_season_std <- ((curr_value_season - season_avg)/(season_rng/2))
          prev_value_season_std <- ((prev_value_season - season_avg)/(season_rng/2))
        }
        else {
          curr_value_season_std <- ((curr_value_season - season_avg))
          prev_value_season_std <- ((prev_value_season - season_avg))
          curr_value_std <- ((curr_value - season_avg))
          prev_value_std <- ((prev_value - season_avg))
        }
        
        immediate_residual <- curr_value_std - prev_value_std
        #immediate_residual <- curr_value - prev_value
        season_residual <- curr_value_season_std - prev_value_season_std
        #season_residual <- curr_value_season - prev_value_season
        
        
        #residual in approximately range [-2,+2]
        residual <- immediate_residual - season_residual
        
        
        #stored for plotting
        cell_residual[[node]][curr_period] <- residual
        cell_val_std[[node]][curr_period] <- curr_value_std
        
        cell_acalc[[node]][curr_period] <- (abs(residual) + 1) ^ 2
        
        
        shortCUSUM <- residual + (prev_season_diff * 0.6)
        longCUSUM <- weighted.mean(c(prev_season_lon_diff, residual), c(0.6,1))
        
        
        cell_season_diff[[node]][curr_period] <- shortCUSUM
        cell_season_lon_diff[[node]][curr_period] <- (abs(shortCUSUM) + 1) ^ 2
        
        #if residual above threshold, OUT : cell changed behavior from past
        if((abs(residual) + 1) ^ 2 > cell_threshold) {
          
          
          #update current state, shift previous for this node in its list
          prev_state[[node]][node,] <- shiftRight(prev_state[[node]][node,] , 1)
          if(residual > 0)
            prev_state[[node]][node,1] <- "OUT-up"
          else
            prev_state[[node]][node,1] <- "OUT-down"
          
          #update cell seasonality
          cell_seasonality[[node]][curr_period] <- weighted.mean(
            c(cell_seasonality[[node]][curr_period],cluster[node,time]),
            c(1,cell_update_rate_OUT))
          
        }
        #else OK
        else {
          
          if((abs(shortCUSUM) + 1) ^ 2 > diff_threshold) {
            #difference is accumulating
            #update current state, shift previous for this node in its list
            prev_state[[node]][node,] <- shiftRight(prev_state[[node]][node,] , 1)
            if(shortCUSUM > 0)
              prev_state[[node]][node,1] <- "TREND-up"
            else
              prev_state[[node]][node,1] <- "TREND-down"
            
            #update cell seasonality
            cell_seasonality[[node]][curr_period] <- weighted.mean(
              c(cell_seasonality[[node]][curr_period],cluster[node,time]),
              c(1,cell_update_rate_OUT))
            
          }
          else{
            #update current state, shift previous for this node in its list
            prev_state[[node]][node,] <- shiftRight(prev_state[[node]][node,] , 1)
            prev_state[[node]][node,1] <- "OK"
            
            #update cell seasonality
            cell_seasonality[[node]][curr_period] <- weighted.mean(
              c(cell_seasonality[[node]][curr_period],cluster[node,time]),
              c(1,cell_update_rate_OK))
          }          
        }
      }    
      
    }
    #END NODE OPERATION
    
    
    
    #Before cluster operation, synch, all nodes finish phase 1
    
    curr_period <- ((time-1) %% period) + 1
    #Propagation phase: each node sends its neighbors within a given range its current status
    #and value
    
    for(node in (1:nrow(cluster))) {
      for(i in (1:nrow(cluster))) {
        if(i != node) {
          prev_state[[node]][i,] <- shiftRight(prev_state[[node]][i,], 1)
          prev_state[[node]][i,1] <- prev_state[[i]][i,1]
        }
        prev_neigh_value[[node]][i,] <- shiftRight(prev_neigh_value[[node]][i,], 1)
        #prev_neigh_value[[node]][i,1] <- cluster[i,time]
        prev_neigh_value[[node]][i,1] <- cell_residual[[i]][curr_period]
      }
    }
    
    for(node in (1:nrow(cluster))) {
      
      
      #START CLUSTER OPERATION
      
      curr_period <- ((time-1) %% period) + 1
      
      
      #VOTING PROCEDURE
      #weighting neighbor to be added
      
      vote <- 0
      
      #set prediction for a single node
#       if(node == watch_node) {
#         if(time > 24)
#           anomalyPredictionOn[node] <- TRUE
#       }
      
      for(nn in (1:neigh_number)) {
     # for(i in (1:nrow(cluster))) {
  
        if(time > 1) {
          debug <- 2
        }
        
        i <- neigh_index[[node]][nn]
        
        ordnum <- as.numeric(ordered(prev_state[[node]][i,1], 
                                     levels = c("OUT-down", "TREND-down", "OK", "TREND-up", "OUT-up")))
        
        #ordnum from -1 to +1
        ordnum <- (ordnum - 3)/2
        
        #count less nodes without anomalies
        if(ordnum > 0)
          ordnum <- ordnum / 2
        
        vote <- vote + neigh_weight[[node]][i]*ordnum
        
        
      }
      
      vote <- vote / neigh_number

      #vote <- sign(vote) * vote^2

      trend_threshold <- 0.6      
      out_threshold <- 0.8
      
      vote_threshold <- out_threshold
      
      
      if(vote < -(vote_threshold)) {
        
        #ANOMALY PROCEDURE
        
        
        if(anomalyPredictionOn[node] == TRUE) { 
          
          #states <- buildStateEntry(node, prev_state, prev_neigh_value, state_mem, nrow(cluster))
          states <- buildStateEntrySeason(node, prev_state, prev_neigh_value, state_mem, nrow(cluster), voteResult = "OK",neigh_weight[[node]])
          
          curr_state <- states[[1]]
          memorized_states <- states[2:state_mem]
          
          anomaly_structure[[node]] <- structure.addNewAnomaly(anomaly_structure[[node]],
                                                               curr_state,
                                                               memorized_states)
          
          
          print(paste("Node ",node,", anomaly time ",time))
          
   
        }
        else {
          
          if(node == watch_node && genAnomaly) {
            
            seasonalValues <- colMeans(prev_neigh_value[[node]])
            seasonalValues <- seasonalValues[state_mem:1]
            
            #seasonalValues <- colMeans(cluster)[(time - state_mem + 1):time]
            
            
            plot(seasonalValues, col="blue", main =paste("Anomaly",time), ylim = c(-2, 2))
            lines(seasonalValues[1:length(seasonalValues) - 1], col = "blue")
            
            diff <- curr_period - state_mem
            if(curr_period - state_mem < 0) {
              cluster_vals <- c(cluster_seasonality[[node]][(length(cluster_seasonality[[node]]) + diff + 1) : length(cluster_seasonality[[node]])],
                                cluster_seasonality[[node]][1:(curr_period)])
            }
            else {
              cluster_vals <- c(cluster_seasonality[[node]][(curr_period - state_mem +1):(curr_period)])
            }
            
            cluster_avg <- mean(cluster_seasonality[[node]], na.rm = T)
            cluster_rng <- (range(cluster_seasonality[[node]])[2] - 
                              range(cluster_seasonality[[node]])[1])
            
            if(cluster_rng != 0) {
              cluster_vals_std <- (cluster_vals - cluster_avg)/cluster_rng
            }
            else{
              cluster_vals_std <- (cluster_vals - cluster_avg)
            }
            lines(cluster_vals_std)
            
            
          }
          print(paste("Node ",node,", anomaly cluster (", paste(neigh_index[[node]], collapse = ",") ,") time ",time))
        }
        
      }
      else {
        
        if(vote < - trend_threshold) {
          
          #detect trend down
          
          print(paste("Node ",node,", trend down time ",time))
          
        }
        
        
        if(anomalyPredictionOn[node] == TRUE) {
          
          #states <- buildStateEntry(node, prev_state, prev_neigh_value, state_mem, nrow(cluster))
          states <- buildStateEntrySeason(node, prev_state, prev_neigh_value, state_mem, nrow(cluster), voteResult = "OK",neigh_weight[[node]])
          
          curr_state <- states[[1]]
          memorized_states <- states[2:state_mem]
          
          newupd <- structure.normalUpdate(anomaly_structure[[node]],
                                           curr_state,
                                           memorized_states)
          
          if(!is.null(newupd))
            anomaly_structure[[node]] <- newupd
          
          anom_hops <- predictAnomaly(anomaly_structure[[node]], 
                                      curr_state,
                                      memorized_states)

            predictAnomaly(anomaly_structure[[node]], 
                           curr_state,
                           memorized_states,
                           look_ahead = 1,
                           genericPrediction = T)
        
          
          anomProbability <- 0
          for(anomaly in anom_hops) {
            anomProbability <- anomProbability + anomaly$hopProb
            if(anomaly$steps > 0)
              print(paste("PredictAnomaly ", anomaly$hopState$id  ," in ",anomaly$steps," with prob ",anomaly$hopProb))
            
          }
          
          if(anomProbability > 0.5) {
            #print("PredictAnomaly within look_ahead = 3")
          }
                    
        }
      }
      
      
      #WEIGHT UPDATE
      for(i in (1:nrow(cluster))) {
        
        neigh_weight[[node]][i] <- (0.4*neigh_weight[[node]][i]) + 
          (0.6*(1 / (1+ (abs(cell_residual[[node]][curr_period] -
                             cell_residual[[i]][curr_period]))^2 )))
        
        #find best neighbors
        #library(FNN)
        #neigh_index[[node]] <- get.knn(neigh_weight[[node]], neigh_number, "kd_tree")[[1]][,]
        
      #  updateProb <- runif(1, 0,1)
        
      #  if(updateProb > 0.7) {
          neigh_index[[node]] <- (head(sort(neigh_weight[[node]], decreasing = T, index.return = T)$ix, neigh_number))
       # }
      }
      
      #Seasonality update (for re-clustering)
      
      #check cluster seasonality
      curr_season_val <- cluster_seasonality[[node]][curr_period]
      
      #if no season value exists, use mean
      if(time < period + 1) {
        
        #update current season value
        cluster_seasonality[[node]][curr_period] <- mean(cluster[,time], na.rm = T)
        
        
      }
      #else correct with formula
      else {
        
        weights <- vector()
        for(i in (1:nrow(cluster))) {
          if(prev_state[[node]][i,1] == "OK") {
            w <- 1
          }
          else {
            w <- 0.3
          }
          weights <- c(weights, w)    
        }
        
        new_season_val <- weighted.mean(cluster[,time], weights, na.rm = T)
        
        cluster_seasonality[[node]][curr_period] <- weighted.mean(
          c(cluster_seasonality[[node]][curr_period],new_season_val),
          c(1,cluster_update_rate))
        
      }  
      
      
      if(predictedAnomaly[[node]] == T) {
        print(paste("Node ",node,"predict next anomaly"))
      }
      
      
      #END CLUSTER OPERATION
    }
    
    stat <- force(paste(prev_state[[watch_node]][watch_node,1], collapse = "-"))
    print(paste("end time",time,", status ", stat))
    
    if(curr_period == 24) {
      #plot(cell_seasonality[[5]], main =paste("cell seasonality, iteration",((time-1) %/% period) +1 ))
      
      plot(cluster[watch_node,c((time - 23) : (time))], col="blue", main =paste("cell seasonality, iteration",((time-1) %/% period) +1 ), ylim = c(1,10))
      lines(cluster[watch_node,c((time - 23) : (time))],col="blue")
      lines(cell_seasonality[[watch_node]],col="black")
      lines(-cell_season_diff[[watch_node]], col="red")
      lines(-cell_residual[[watch_node]], col = "green")
      points(cell_acalc[[watch_node]], col = "red")
      points(cell_season_lon_diff[[watch_node]], col = "plum3")
    }
  }
  
  plot(cluster_seasonality[[watch_node]], main = "Final cluster seasonality")
  
  print("end")
  
  
}

#shift cells of a vector one position to the right, with wrapping
shiftRight <- function(vec, num = 1) {
  
  vec <- c(vec[(length(vec) - num + 1):length(vec)],vec[1:(length(vec)-num)])
  return(vec)
  
}

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


#returns a list containing next p steps with attached probability according to the cell historical data given
#MST = Markov State Transition
cellPrevisionMarkov <- function(cellFeatures, cellMST) {

  
  
  
  
}

