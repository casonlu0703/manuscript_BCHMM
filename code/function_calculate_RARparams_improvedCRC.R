# Given a candidate start hour of a day, find
# (1) the center resting clock (CRC),
# (2) the cumulative state probability centered around the CRC.
center_clock_index <- function(
    starthr        # start hour of a 24-hour day, values = 0, 4, 8, 12, 16, 20, i.e., 0 am to 8 pm, every 4 hours
    , state        # specify which hidden activity state, values = 1, 2, 3, i.e., inactive, moderately active, and highly active
    , one.day.prob # state probabilities for multiple 24-hour days
){
  # rm(list = c("starthr", "state", "one.day.prob"))
  
  sf <- 1/12  #sampling frequency in hour unit, i.e. here sf = 5min/60min = 1/12 
  clock.time <- seq(0, 24*2, sf)
  
  # 24 hours from (starthr) to (starthr + 24)
  # 1 = inactive state probability 
  # 2 = moderately active
  # 3 = highly active
  # NOTICE: length should always represents one 24 hour
  # for example, our data are collected in 5 mins. The length is therefore 288 = 24 * 60/5
  p <- one.day.prob[
    (which(clock.time == starthr)):(which(clock.time == starthr + 24) - 1)
    , state
  ]
  
  # for inactivity: amount of rest, i.e., duration of rest per day (hours) 
  # for highly active: amount of activity, i.e., duration of being highly active (hours)
  # rest_amount is invariant to the definition of start hour of a day
  rest_amount <- 24*sum(p)/length(p)
  
  # absolute index to compute the gravity center of p1
  index <- seq(0, 24 - 1/12, 1/12)  
  
  # center of rest which corresponds to the gravity center of p1
  # Note: only valid with suitable hour_day_start where the major sleep/rest is of one- phase. 
  center_rest <- sum(p * index / sum(p)) 
  
  # center_rest is the position in index not the clock time. 
  # One can easily convert it to clock time via: 
  index_clocktime <- seq(starthr, starthr + 24 - sf, sf) 
  
  #clocktime[index_clocktime > 24] <- (index_clocktime - 24)[index_clocktime > 24] 
  center_rest_clock <- index_clocktime[which.min(abs(index - center_rest))] 
  
  # cumulative probability centered at the center resting clock
  resting_clock_start <- clock.time[which.min(abs(clock.time - (center_rest_clock - rest_amount/2) ))]
  resting_clock_end   <- clock.time[which.min(abs(clock.time - (center_rest_clock + rest_amount/2) ))]
  
  cumP <- sum(
    one.day.prob[
      (which(clock.time == resting_clock_start)):(which(clock.time == resting_clock_end))
      , state
    ])
  
  # the clock time of the highest state probability 
  MaxPP <- index_clocktime[which(p == max(p))]
  
  #worst clock : complete lack of circadian rhythm where the probability of rest is constant and equal to rest_amount/24
  #worst_p1 <- rep(mean(p1), 24/sf) 
  
  # perfect clock: rest period with no interruptions 
  find_perfect_p <- function(center_rest, rest_amount, index, sf){ 
    perfect_p <- rep(0, 24/sf) 
    t1 <- center_rest - rest_amount/2 
    t2 <- center_rest + rest_amount/2 
    perfect_t2 <- which.min(abs(index - center_rest)) 
    perfect_t1 <- which.min(abs(index - center_rest + rest_amount/2)) 
    perfect_t3 <- which.min(abs(index - center_rest - rest_amount/2)) 
    perfect_p[perfect_t1:perfect_t3] <- 1 
    if(t2 > max(index)){
      offset <- round((t2-max(index))/sf)
      perfect_p[1:offset]<-1
    } 
    if(t1 < min(index)){
      offset <- round((min(index)-t1)/sf)
      perfect_p[-offset+L, L]<-1
    } 
    return(perfect_p) 
  } 
  
  perfect_p <- find_perfect_p(center_rest, rest_amount, index, sf) 
  
  #rhythm index RI, see in eq(4) of (Huang, 2018). 
  RI <- (sum(p[which(perfect_p > 0)])*sf/rest_amount - rest_amount/24)*24/(24 - rest_amount) 
  
  return(c(starthr, cumP, center_rest_clock, rest_amount, MaxPP, RI))
}


circadian.param <- function(
    one.day.prob # the state probabilities
    ){
  
  # Calculate RAR parameters for IA (state 1) -----
  # candidate start hour of a 24-hour day
  starthr_list <- seq(0, 20, by = 4) # define the start of a one-day state probabilities

  CRC_list <- as.data.frame(t(as.data.frame(
    lapply(
      starthr_list
      , FUN = center_clock_index
      , state = 1, one.day.prob = one.day.prob
    )
  )))
  
  rownames(CRC_list) <- NULL
  colnames(CRC_list) <- c("StartHr", "cumulative_StateProb", "CRC", "RA", "MIPP", "RI")
  
  # convert to 24 hours
  CRC_list$CRC[CRC_list$CRC > 24] <- CRC_list$CRC[CRC_list$CRC > 24] - 24
  CRC_list$MIPP[CRC_list$MIPP > 24] <- CRC_list$MIPP[CRC_list$MIPP > 24] - 24
  
  CRC_list_1 <- CRC_list[which.max(CRC_list$cumulative_StateProb), ]
  
  # Calculate RAR parameters for HA (state 3) -----
  # candidate start hour of a 24-hour day
  starthr_list <- seq(0, 20, by = 4) # define the start of a one-day state probabilities
  
  CAC_list <- as.data.frame(t(as.data.frame(
    lapply(
      starthr_list
      , FUN = center_clock_index
      , state = 3, one.day.prob = one.day.prob
    )
  )))
  
  rownames(CAC_list) <- NULL
  colnames(CAC_list) <- c("StartHr", "cumulative_StateProb", "CAC", "AA", "MAPP", "RI_active")
  
  # convert to 24 hours
  CAC_list$CAC[CAC_list$CAC > 24] <- CAC_list$CAC[CAC_list$CAC > 24] - 24
  CAC_list$MAPP[CAC_list$MAPP > 24] <- CAC_list$MAPP[CAC_list$MAPP > 24] - 24
  
  CAC_list_1 <- CAC_list[which.max(CAC_list$cumulative_StateProb), ]
  
  return(
    list(
      S1StartHr   = CRC_list_1$StartHr
      , RI        = CRC_list_1$RI
      , RA        = CRC_list_1$RA
      , CRC       = CRC_list_1$CRC
      , MIPP      = CRC_list_1$MIPP
      , S3StartHr = CAC_list_1$StartHr
      , RI_active = CAC_list_1$RI_active
      , AA        = CAC_list_1$AA
      , CAC       = CAC_list_1$CAC
      , MAPP      = CAC_list_1$MAPP
      )
    )

}
