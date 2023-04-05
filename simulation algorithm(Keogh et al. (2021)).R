############################################################
#Simulation algorithm proposed by Keogh et al. (2021) #
############################################################
library(locfit)#for inverse logit function
library(glue)

#Simulation function
"""
Input:
TT:the number of time points to follow up, a numerical value (default = 4, means 5 visits:k=0,1,2,3,4)
alpha: a vector of paramters for conditional hazard model, of length of 4
n: the number of subjects to simulate

Output:
Subject_id: the id of subject
Time: observed time point(visit)
U: unmeasured covariate (individual frailty) values acorss time
L: measured covariate values (biomarker) across time
A: binary treatment indicator across time
Alag: lag values of binary treatment indicator across time
Hazard: conditional hazard estimated by a given conditional hazard model across time
Time_obs: observed event time
D_obs: event indicator of individual(1: event happened for this individual; 0: event didn't happen for this individual)
event: event indicator of individual at certain visit(1: event happend at this time point; 0:event didn't happen at this time point)
time_stop: the time to stop the follow-up at certain visit time
"""
sim_algorithm_Keogh = function(TT,alpha,n){
  output_df = data.frame()
  for (i in 1:n){
    U = vector(mode = 'numeric',length = 1)
    L = vector(mode = 'numeric',length = 1)
    A = vector(mode = 'numeric',length = 1)
    hazard = vector(mode = 'numeric',length = 1)
    T_obs = NA
    event = vector(mode = 'numeric', length = 1)
    #Initialization
    U = rnorm(1,mean = 0,sd = 0.1)
    L[1] = rnorm(1,mean = U[1],sd = 1)
    A[1] = rbinom(n = 1,size = 1,prob = expit(-2+0.5*L[1]))
    #compute hazard from a given conditional hazard in the period 0<t<1, to determine the event time T_obs
    hazard[1] = alpha[1]+alpha[2]*A[1]+alpha[3]*L[1]+alpha[4]*U[1]
    V = runif(n = 1,min = 0,max = 1)
    delta_T = -log(V)/hazard[1]
    T_obs = ifelse(is.na(T_obs) & delta_T<1 & hazard[1]>0,delta_T, T_obs)#using additive hazard model is likely to generate hazard with negative values with tiny possibility
    #if no event occurred, continue sampling the data
    for (t in 2:(TT+1)){
      if (is.na(T_obs)){
        L[t] = rnorm(1, mean = 0.8*L[t-1]-A[t-1]+0.1*t+U[1], sd = 1)
        A[t] = rbinom(n = 1,size = 1, prob = expit(-2+0.5*L[t]+A[t-1]))
        #compute hazard from a given conditional hazard in the period s<=t<s+1, to determine the event time T_obs
        hazard[t] = alpha[1]+alpha[2]*A[t]+alpha[3]*L[t]+alpha[4]*U[1]
        V = runif(n = 1,min = 0,max = 1)
        delta_T = -log(V)/hazard[t]
        T_obs = ifelse(is.na(T_obs) & delta_T<1 & hazard[t]>0,t-1+delta_T, T_obs)#using such addtive hazard is likely to generate hazard with negative values with tiny possibility; t-1 is exactly "s" in the pseudo-code
      }else{
        break
      }
    }
    #if no event occurred in the period 0<t<TT+1, subject i is administratively censored at time t = TT + 1
    #ex:under TT = 4 (5 visits), if no event happened within 0<t<5, then censored at time 5
    T_obs=ifelse(is.na(T_obs),TT+1,T_obs)
    #longitudinal data for subject i:
    subject_id = rep(i,length(L))
    Time = seq(0,length(L)-1)
    Time_obs = rep(T_obs,length(L))
    D_ind = ifelse(T_obs<TT+1, 1, 0)
    D_obs = rep(D_ind,length(L))#event indicator of individual i 
    time_stop=Time+1
    time_stop=ifelse(time_stop>Time_obs,Time_obs,time_stop) #time_stop: used to find the time point of event 
    event=ifelse(time_stop==Time_obs & D_obs==1,1,0)#event indicator of individual i at certain time point
    df = data.frame(Subject_id = subject_id,
                    Time = Time,
                    Time_obs = Time_obs,
                    D_obs = D_obs,
                    event = event,
                    time_stop = time_stop,
                    Hazard = hazard,
                    U = U,
                    L = L,
                    A = A
    )
    #add lag values of treatment binary indicator A to the dataframe
    for (j in 1:TT){
      end_idx = length(A)-j
      if (end_idx>0){
        df$Alag = c(rep(0,j),A[1:end_idx])
      }else{
        df$Alga = rep(0,length(A))
      }
      colnames(df)[ncol(df)] = glue("Alag{j}")
    }
    output_df = rbind(output_df, df)
  }
  return (output_df)
}

#ex:simulating longitudinal data for n = 100 subjects 
alpha = c(0.7,-0.2,0.05,0.05)
TT = 4 #indicates 5 visits: k = 0,1,2,3,4
n = 100 #100 individuals
data = sim_algorithm_Keogh(TT , alpha, n)
View(data)
    
  