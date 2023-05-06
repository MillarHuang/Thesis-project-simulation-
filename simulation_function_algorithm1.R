############################################################
#Simulation algorithm proposed by Havercroft et al. (2012) #
############################################################

#Simulation function
# """
# Input:
# TT:the number of time points to follow up, a numerical value
# k: the check-up time, a numerical value
# gamma: a vector of paramters for MSM hazard, of length of 4
# theta: a vector of paramters for logistic regression in treatment assignment, of length of 3
# n: the number of subjects to simulate
# 
# Output:
# Subject_id: the id of subject
# Time: observed time points
# U: unmeasured covariate values acorss time
# L: measured confounding(CD4 count) across time
# A: binary treatment indicator across time
# Hazard: hazard estimated by a given MSM across time
# Y_next_time_point: binary failure indicator for next time point(t+1)
# Time_to_treament: treament initialiaztion time (if it equals to NA, it means subject hasn't recevied any treatment)
# 
# """
sim_algorithm = function(TT,k,gamma,theta,n,tau){
  output_df = data.frame()
  for (i in 1:n){
    U = vector(mode = 'numeric',length = 1)
    L = vector(mode = 'numeric',length = 1)
    A = vector(mode = 'numeric',length = 1)
    Y = vector(mode = 'numeric',length = 1)
    hazard = vector(mode = 'numeric',length = 1)
    T_init = TT+1 #initialize as T+1: if final T_init is smaller or equal to T, it indicates that subject received treatment during the T time points, otherwise didn't receive any treatment
    #Initialization
    U[1] = runif(1,min=0,max=1)
    #simulate baseline CD4 counts(L_0)
    error = rnorm(1,0,20)
    L[1] = qgamma(U[1],shape = 3,scale = 154) + error
    #draw treatment decision(A_0)
    ##############################################
    #*Insert violation of positivity assumption: 
    #*when L_t <= tau & A_t-1 = 0 (A_-1 = 0), the subject is forced to start the treatment
    ################################################
    if (L[1]<= tau){
      A[1] = 1
    }else{
      prob = expit(theta[1] + theta[3]*(L[1]-500))
      A[1] = rbinom(n = 1,size = 1,prob = prob)
    }
    #if treatment starts at time 0, then set T_init as 0
    if (A[1] == 1){
      T_init = 0
    }
    #calculate the hazard from a given MSM, to generate the binary failure indicator
    hazard[1] = expit(gamma[1]+gamma[3]*A[1])
    if(hazard[1] >= U[1]){
      Y[1] = 1
    }else{
      Y[1] = 0
    }
    # 2:T+1 represents time point 1 to T
    for (t in 2:(TT+1)){
      if(Y[t-1] == 0){
        #update Ut
        noise = rnorm(1,0,0.05)
        U[t] = min(1,max(0,U[t-1]+noise))
        #not at a check-up time, remains previous values
        if ((t-1) %% k != 0){
          L[t] = L[t-1]
          A[t] = A[t-1]
        }else{
          #at a check_up time
          #if treatment starts at the last check-up, CD4 count will have a shift(+150)
          error = rnorm(1,100*(U[t]-2),50)
          if((t-1)%/%k == 1){
            #if at the second check-up(t0 is the first check-up), t-k-1 is 0
            L[t] = max(0,L[t-1]+150*A[t-k]+error)
          }else{
            L[t] = max(0,L[t-1]+150*A[t-k]*(1-A[t-k-1])+error)
          }
          #Once starts treatment, it will remain on treatment until failure or end of study
          ##############################################
          #*Insert violation of positivity assumption: 
          #*when L_t <= tau & A_t-1 = 0, the subject is forced to start the treatment
          ################################################
          if(A[t-1] == 0){
            if(L[t] <= tau){
              A[t] = 1
            }else{
              prob = expit(theta[1] + theta[2]*(t-1) + theta[3]*(L[t]-500))
              A[t] = rbinom(n = 1,size = 1,prob = prob)
            }
          }else{
            A[t] = 1
          }
          #if treatment starts at this time point, set T_init as this time point
          if(A[t] == 1 & A[t-k] == 0){
            T_init = t-1
          }
        }
        #compute the hazard to generate the binary failure indicator
        T_star = min(t-1,T_init)#if treatment hasn't start yet, T^* in the formula is equal to the present time
        hazard[t] = expit(gamma[1]+gamma[2]*((1-A[t])*(t-1)+A[t]*T_star)+gamma[3]*A[t]+gamma[4]*A[t]*((t-1)-T_star))
        if( (1-prod(1-hazard)) >= U[1]){
          Y[t] = 1
          
        }else{
          Y[t] = 0
        }
      }else{
        break
      }
    }
    #longitudinal data for subject i
    subject_id = rep(i,length(U))
    Time = seq(0,length(U)-1)
    Time_to_treatment = ifelse(T_init <= TT, rep(T_init,length(U)),NA)
    df = data.frame(Subject_id = subject_id,
                    Time = Time,
                    U = U,
                    L = L,
                    A = A,
                    Hazard = hazard,
                    Y_next_time_point = Y,
                    Time_to_treatment = Time_to_treatment
    )
    #add lag values of A (lag==1)
    lag=1
    for (j in 1:lag){
      end_idx = length(A)-j
      if (end_idx>0){
        df$Alag = c(rep(0,j),A[1:end_idx])
      }else{
        df$Alag = rep(0,length(A))
      }
      colnames(df)[ncol(df)] = glue("Alag{j}")
    }
    output_df = rbind(output_df, df)
  }
  return(output_df)
}


#Simulation algorithm default parameters setting
TT=40
k=5
gamma = c(-3, 0.05, -1.5, 0.1)
theta = c(-0.405, 0.0205, -0.00405)

