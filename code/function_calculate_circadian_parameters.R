circadian.param <- function(one.day.prob){
  #Most healthy people sleep/rest at night so generally we choose 12pm-11:59am to present a day (with one-phase sleep/rest): 
  hour_day_start<-12 
  # The selection of hour_day_start will affect center rest time and Rhythm Index (RI). I suggest one should choose hour_day_start which leads to highest RI and by experience, hour_day_start =12 is suitable for most healthy subjects. 
  #For late-type person, one should select another suitable time segment, i.e. 0am- 23:59pm 
  library('lubridate') 
  # this function is used to find an appropriate time index of a day start 
  find_day_start<-function(hour_day_start,time){ 
    hour_time<-hour(time)
    min_time<-minute(time) 
    A<-which(hour_time==hour_day_start) 
    B<-which(min_time<60*sf) 
    one_day_start<-A[min(which(A%in%B, arr.ind = TRUE))] 
    return(one_day_start) 
  } 
  
  #one_day_prob <- as.data.frame(result1[[1]][[2]][289:(289+288-1),])
  #one_day_prob <- sp0[289:(289+288-1),1:3]
  one_day_prob <- one.day.prob
    
  p1<-one_day_prob[,1] #inactive state probability 
  p2<-one_day_prob[,2] #moderately active state probability 
  p3<-one_day_prob[,3] #highly active state probability 
  
  ########################################################## 
  #Following circadian parameters are calculated based on one_day_prob
  #see details in Circadian Parameters Section of the main manuscript 
  ########################################################## 
  #amount of rest: duration of rest per day (hours) 
  rest_amount<-24*sum(p1)/nrow(one_day_prob) 
  index<-seq(0,24-1/12,1/12) #absolute index to compute the gravity centre of p1 #center of rest which corresponds to the gravity centre of p1
  #Note: only valid with suitable hour_day_start where the major sleep/rest is of one- phase. 
  center_rest<-sum(p1*index/sum(p1)) 
  #center_rest is the position in index not the clock time. #One can easliy convert it to clock time via: 
  sf<-1/12 #sampling frequency in hour unit, i.e. here sf=5min/60min=1/12 
  index_clocktime<-seq(hour_day_start,hour_day_start+24-sf,sf) 
  clocktime<-index_clocktime 
  clocktime[index_clocktime>24]<-(index_clocktime-24)[index_clocktime>24] 
  center_rest_clock<-clocktime[which.min(abs(index-center_rest))] 
  #worst clock : complete lack of circadian rhythm where the probability of rest is constant and equal to rest_amount/24
  worst_p1<-rep(mean(p1),24/sf) 
  
  #perfect clock: rest period with no interruptions 
  find_perfect_p1<-function(center_rest,rest_amount,index,sf){ 
    perfect_p1<-rep(0,24/sf) 
    t1<-center_rest-rest_amount/2 
    t2<-center_rest+rest_amount/2 
    perfect_t2<-which.min(abs(index-center_rest)) 
    perfect_t1<-which.min(abs(index-center_rest+rest_amount/2)) 
    perfect_t3<-which.min(abs(index-center_rest-rest_amount/2)) 
    perfect_p1[perfect_t1:perfect_t3]<-1 
    if(t2>max(index)){offset<-round((t2-max(index))/sf);perfect_p1[1:offset]<-1} 
    if(t1<min(index)){offset<-round((min(index)-t1)/sf);perfect_p1[-offset+L,L]<-1} 
    return(perfect_p1) 
  } 
  perfect_p1<-find_perfect_p1(center_rest,rest_amount,index,sf) 
  #rhythm index RI, see in eq(4) of the main manuscript. 
  RI<-(sum(p1[which(perfect_p1>0)])*sf/rest_amount-rest_amount/24)*24/(24- rest_amount) 
  
  return(list(RI=RI, RA=rest_amount, CRC=center_rest_clock))
  ########################################################## 
  #the function figure_day_profile generates day profile plots like Figure 7 in the manuscript. ########################################################## 
  #We choose 12pm-11:59am to present one day so that use the following setting in the plot. For other hour_day_start, the following parameters need to be adjusted. 
  # index_position<-seq(0,24,2)
  # index_lable<-c('12','14','16','18','20','22', 
  #                '0/24','2','4','6','8','10','12') #corresponding clock time #Note: one need to change index_lable if the day segment changes. 
  # source('../../../HMMS-master/figure_day_profile.R') 
  # par(fig=c(0,1,0,1), new=F) 
  # figure_day_profile(id=0,one_day_prob,index,index_position,index_lable) 
  # figure_day_profile(id=1,as.data.frame(result1[[1]][[2]][289:(289+288-1),]),index,index_position,index_lable) 
  
}
