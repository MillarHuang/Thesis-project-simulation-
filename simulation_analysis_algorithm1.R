#set seed
set.seed(10)

#the MSM parameters to estimate
gamma1=vector(mode = 'numeric',length = 100)
gamma2=vector(mode = 'numeric',length = 100)
gamma3=vector(mode = 'numeric',length = 100)

B = 100
#B = 100 replications each of a study with n subjects
for(i in 1:B){
  #simulate longitudinal data
  sim_data = sim_algorithm(TT=TT,k=k,gamma=gamma,theta=theta,n=n,tau=tau)
  #------------------
  #Stabilized weights
  #------------------
  #use the subset that before initializing treatment(I(A_t−1 = 0)) & the data at check-up time: to fit P(At|At-1,L) and P(At|At-1)
  sim_data_bt = sim_data[sim_data$Alag1==0 & sim_data$Time%%k==0,]
  #denominator
  wt.mod=glm(A~Time+L,family="binomial",data=sim_data_bt)
  pred.wt=predict(wt.mod,type = "response")
  sim_data_bt$wt=ifelse(sim_data_bt$A== 1,pred.wt,1-pred.wt)
  sim_data_bt$wt.cum=ave(sim_data_bt$wt,sim_data_bt$Subject_id,FUN=cumprod)
  #numerator
  wt.mod.num=glm(A~Time,family="binomial",data=sim_data_bt)
  pred.wt.num=predict(wt.mod.num,type = "response")
  sim_data_bt$wt.num=ifelse(sim_data_bt$A==1,pred.wt.num,1-pred.wt.num)
  sim_data_bt$wt.cum.num=ave(sim_data_bt$wt.num,sim_data_bt$Subject_id,FUN=cumprod)
  #stabilized weight under time-dependent treatment (cumulative)
  sim_data_bt$ipw.s=sim_data_bt$wt.cum.num/sim_data_bt$wt.cum
  #add the stabilized weight to the dataset
  sim_data = left_join(sim_data,data.frame(Subject_id = sim_data_bt$Subject_id,
                                           Time = sim_data_bt$Time,
                                           ipw.s = sim_data_bt$ipw.s),
                       by=c('Subject_id',"Time"))
  #Inpute the missing values in ipw.s by its previous values, since:
  #between the check-up time, since the A and L do not update, so P(At(=At-1)|At-1) =1, P(At(=At-1)|At-1 ,L)) = 1; 
  #after initializing treatment, P(At = 1|At-1 = 1) = 1, P(At =1 |At-1 = 1,L)) = 1
  sim_data$ipw.s = na.locf(sim_data$ipw.s)
  
  #Set up d1 and d3
  d1 = map2(sim_data$Time,sim_data$Time_to_treatment,~min(.x,.y,na.rm=TRUE))
  d1 = unlist(d1, use.names=FALSE)
  sim_data$d1 = d1
  d3 = map2(sim_data$Time-sim_data$Time_to_treatment,0,~max(.x,.y,na.rm=TRUE))
  d3 = unlist(d3, use.names = FALSE)
  sim_data$d3 = d3
  
  #Fit MSM using IPTW
  msm=glm(Y_next_time_point~d1+A+d3,family = 'binomial',weights = ipw.s,data = sim_data)
  gamma1[i] = coef(msm)['d1']
  gamma2[i] = coef(msm)['A']
  gamma3[i] = coef(msm)['d3']
}
#########################################
#Simulation results
# mean(gamma1)#True:0.05,we get 0.051
# sd(gamma1)#0.011
# mean(gamma2)#True:-1.50,we get -1.512
# sd(gamma2)#0.154
# mean(gamma3)#True:0.10,we get 0.101
# sd(gamma3)#0.012
