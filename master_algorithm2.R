#--------------------------
#packages
#--------------------------

library(survival)
library(timereg)#for additive hazard model :aalen()
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(glue)
library(patchwork)

#--------------------------
# functions.R
#--------------------------
# - Loads the functions needed.

source("simulation_function_algorithm2.R")


#--------------------------
# simulation_analysis_True_algorithm2.R
#--------------------------
# - Obtains (estimated) true values of the cumulative coefficients of the correctly-specified MSM, 
# and of the causal estimands of interest by fitting the MSM to the large trial data generated in dat_sim_TRUE_VALUES_1.R. 
# **Note:it will take multiple hours to run (for B = 1000 repetitions)

source("simulation_analysis_True_algorithm2.R")


#generate results under different combination of sample size(n) and positivity assumption violation threshold(tau)
n_ = c(100,200,300,500,1000)
tau_ = c(-0.3,0.3,0.6,1,1.5) # approximate 0.5,0.6,0.7,0.8,0.9 quantile of L
#calculate the IPTW results and generate the corresponding table, then save data as csv file
for (n in n_){
  for (tau in tau_){
    #obtain IPTW estimates
    source("simulation_analysis_algorithm2.R")
    #generate two table2 and saved as csv file: 1. cumulative coefficients for each parameters 2.Survival probabilities for the treatment regimes ‘never treated’ and ‘always treated
    source("simulation_results_algorithm2.R")
  }
}

#import and combine the results
results_table1 = data.frame()
results_table2 = data.frame()
for (n in n_){
  for(tau in tau_){
    results_long1 = read_csv(glue("results_data/violation_result_table1_algorithm2_{tau}_{n}.csv"))
    results_table1 = rbind(results_table1,results_long1)
    results_long2 = read_csv(glue("results_data/violation_result_table2_algorithm2_{tau}_{n}.csv"))
    results_table2 = rbind(results_table2,results_long2)
  }
}
write.table(results_table1,"results_data/violation_result_table1_algorithm2_all.csv",sep=",",row.names = F)
write.table(results_table2,"results_data/violation_result_table2_algorithm2_all.csv",sep=",",row.names = F)



#import results
results_table1 = read_csv("results_data/violation_result_table1_algorithm2_all.csv")
results_table2 = read_csv("results_data/violation_result_table2_algorithm2_all.csv")

#######################################################################################
#table 1
#######################################################################################
#Split the strings in columns Mean_estimate__Empirical_SE and Bias__MonteCarlo_SE of dataset
ff_split = function(x){
  output = str_split(x," ",simplify= TRUE)
  output[2] = gsub("[\\(\\)]", "", regmatches(output[2], gregexpr("\\(.*?\\)", output[2]))[[1]])
  return(as.numeric(output))
}
mean_se = matrix(unlist(map(results_table1$Mean_estimate__Empirical_SE,ff_split)),ncol=2,byrow = TRUE)
bias_se = matrix(unlist(map(results_table1$Bias__MonteCarlo_SE,ff_split)),ncol=2,byrow = TRUE)
results_table1$Mean_estimate = mean_se[,1]
results_table1$Empirical_SE = mean_se[,2]
results_table1$Bias = bias_se[,1]
results_table1$MonteCarlo_SE = bias_se[,2]
results_table1$Mean_estimate__Empirical_SE=NULL
results_table1$Bias__MonteCarlo_SE = NULL

#Combine the True values in column Mean_estimate, while set its tau as 'True' factor level
cum_coef = c("alpha_0(t)","alpha_A,0(t)", "alpha_A,1(t)", "alpha_A,2(t)" ,"alpha_A,3(t)", "alpha_A,4(t)")
time = length(unique(results_table1$Time))
for(i in cum_coef){
  for (nn in n_){
    tt = results_table1%>%filter(Cumulative_coefficient == i,n == nn)%>%head(time)
    tt$tau = 'True'
    tt$Mean_estimate = tt$True_value
    tt$Empirical_SE = 0
    tt$Bias = 0
    tt$MonteCarlo_SE = 0
    results_table1= rbind(results_table1,tt)
  }
}

#Plotting
#####################
#IPTW estimator mean#
#####################
#plot function for the data according to a specific cumulative_coefficient
make_plot = function(data, cumulative_coefficient){
  data%>%
  ggplot(aes(x = Time,y = Mean_estimate,color = factor(tau)))+
  geom_point(aes(shape = factor(tau)))+
  geom_line(aes(linetype = factor(tau)))+
  scale_x_continuous(
    breaks = seq(1,5,1),
    labels = seq(1,5,1),
    limits = c(1,5),
    name = 'Time'
  )+
  scale_y_continuous(
    name = glue("{cumulative_coefficient}")
  )+
  scale_color_manual(
      values = c("black","#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
  guides(color=guide_legend(title="Tau"),
         shape = guide_legend(title="Tau"),
         linetype = guide_legend(title="Tau"))
}
#The list of all combined plots(p_m) for different sample size
p_m = list()
#create combined plots for each sample size (n)
for(i in seq(1,length(n_))){
  #generate plot for each cumulative coefficient under a certain sample size n
  plots_ = 
    results_table1 %>%
    #reorder the factor level
    mutate(tau = fct_relevel(tau,'True','-0.3','0.3','0.6','1','1.5'))%>%
    nest(data = -c(n,Cumulative_coefficient)) %>%
    filter(n==n_[i])%>%
    mutate(
      plots = map2(data, Cumulative_coefficient, make_plot)
    ) %>%
    pull(plots) 
  #combine the plots of different cumulative coefficient together for each certain sample size n 
  p_m[[i]] = 
    (plots_[[1]] | plots_[[2]]) / (plots_[[3]] | plots_[[4]])/(plots_[[5]] | plots_[[6]])  +
    plot_annotation(
      subtitle = glue("n = {n_[i]}"),
      title = "IPTW estimator mean under violation of positivity assumption",
    ) &
    theme_classic(8)+
    theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
          plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
}

#save the plots
for(i in seq(1,length(p_m))){
  ggsave(p_m[[i]],filename=glue("results_plot/figure_mean_algorithm2_{n_[i]}.png"),width = 7, height = 6) 
}

##############
#Empirical SE#
##############
make_plot_se = function(data, cumulative_coefficient){
  data%>%
    ggplot(aes(x = Time,y = Empirical_SE,color = factor(tau)))+
    geom_point(aes(shape = factor(tau)))+
    geom_line(aes(linetype = factor(tau)))+
    scale_x_continuous(
      breaks = seq(1,5,1),
      labels = seq(1,5,1),
      limits = c(1,5),
      name = 'Time'
    )+
    scale_y_continuous(
      name = glue("{cumulative_coefficient} SE")
    )+
    scale_color_manual(
      values = c("#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
    guides(color=guide_legend(title="Tau"),
           shape = guide_legend(title="Tau"),
           linetype = guide_legend(title="Tau"))
}
#The list of all combined plots(p_m) for different sample size
p_m_se = list()
#create combined plots for each sample size (n)
for(i in seq(1,length(n_))){
  #generate plot for each cumulative coefficient under a certain sample size n
  plots_ = 
    results_table1 %>%
    filter(tau!='True')%>%
    #reorder the factor level
    mutate(tau = fct_relevel(tau,'-0.3','0.3','0.6','1','1.5'))%>%
    nest(data = -c(n,Cumulative_coefficient)) %>%
    filter(n==n_[i])%>%
    mutate(
      plots = map2(data, Cumulative_coefficient, make_plot_se)
    ) %>%
    pull(plots) 
  #combine the plots of different cumulative coefficient together for each certain sample size n 
  p_m_se[[i]] = 
    (plots_[[1]] | plots_[[2]]) / (plots_[[3]] | plots_[[4]])/(plots_[[5]] | plots_[[6]])  +
    plot_annotation(
      subtitle = glue("n = {n_[i]}"),
      title = "IPTW estimator empirical standard error under violation of positivity assumption",
    ) &
    theme_classic(8)+
    theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
          plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
}
#save the plots
for(i in seq(1,length(p_m_se))){
  ggsave(p_m_se[[i]],filename=glue("results_plot/figure_se_algorithm2_{n_[i]}.png"),width = 7, height = 6) 
}



##############################################################################
#table 2
##############################################################################
#Split the strings in columns Mean_estimate__Empirical_SE and Bias__MonteCarlo_SE of dataset
mean_se = matrix(unlist(map(results_table2$Mean_estimate__Empirical_SE,ff_split)),ncol=2,byrow = TRUE)
bias_se = matrix(unlist(map(results_table2$Bias__MonteCarlo_SE,ff_split)),ncol=2,byrow = TRUE)
results_table2$Mean_estimate = mean_se[,1]
results_table2$Empirical_SE = mean_se[,2]
results_table2$Bias = bias_se[,1]
results_table2$MonteCarlo_SE = bias_se[,2]
results_table2$Mean_estimate__Empirical_SE=NULL
results_table2$Bias__MonteCarlo_SE = NULL

#Combine the True values in column Mean_estimate, while set its tau as 'True' factor level
treat_regime = c("Never treated", "Always treated")
for(i in treat_regime){
  for (nn in n_){
    tt = results_table2%>%filter(Treatment_regime == i,n == nn)%>%head(time)
    tt$tau = 'True'
    tt$Mean_estimate = tt$True_value
    tt$Empirical_SE = 0
    tt$Bias = 0
    tt$MonteCarlo_SE = 0
    results_table2= rbind(results_table2,tt)
  }
}
#Plotting
####################################
#survival probability mean estimate#
####################################
#plot function for the data according to a specific sample size(n): including curve of survival probabilities of Always treated & Never treated
make_plot2 = function(data){
  data%>%
    ggplot(aes(x = Time,y = Mean_estimate,color = factor(tau)))+
    geom_point(aes(shape = factor(tau)))+
    geom_line(aes(linetype = factor(Treatment_regime)))+
    scale_x_continuous(
      breaks = seq(1,5,1),
      labels = seq(1,5,1),
      limits = c(1,5),
      name = 'Time'
    )+
    scale_y_continuous(
      limits = c(0,0.7),
      breaks = seq(0,0.7,0.1),
      name = glue("Survival probabilities")
    )+
    scale_color_manual(
      values = c("black","#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
    guides(color=guide_legend(title="Tau"),
           shape = guide_legend(title="Tau"),
           linetype = guide_legend(title="Treatment regime"))
}

# #The list of all plots(p_m2) for different sample size
# p_m2 = list()
# #create plots for each sample size (n)
# for(i in seq(1,length(n_))){
#   #generate plot for each sample size n
#   plots_ = 
#     results_table2 %>%
#     #reorder the factor level
#     mutate(tau = fct_relevel(tau,'True','-0.3','0.3','0.6','1','1.5'))%>%
#     nest(data = -n) %>%
#     filter(n==n_[i])%>%
#     mutate(
#       plots = map(data, make_plot2)
#     ) %>%
#     pull(plots) 
#   #entitle the plot
#   p_m2[[i]] = 
#     plots_[[1]] +
#     plot_annotation(
#       subtitle = glue("n = {n_[i]}"),
#       title = "IPTW estimator mean of survival probabilities under violation of positivity assumption",
#     ) &
#     theme_classic(8)+
#     theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
#           plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
# }
# #save the plots
# for(i in seq(1,length(p_m2))){
#   ggsave(p_m2[[i]],filename=glue("results_plot/figure_survival_probabilities_algorithm2_{n_[i]}.png"),width = 7, height = 6) 
# }

#Or combine all plots in one plot and save as png file
p_m2 = list()
for(i in seq(1,length(n_))){
  #generate plot for each sample size n
  plots_ = 
    results_table2 %>%
    #reorder the factor level
    mutate(tau = fct_relevel(tau,'True','-0.3','0.3','0.6','1','1.5'))%>%
    nest(data = -n) %>%
    filter(n==n_[i])%>%
    mutate(
      plots = map(data, make_plot2)
    ) %>%
    pull(plots) 
  #add the plot
  p_m2[[i]] = plots_[[1]]
  }
#combine the plots in one plot    
p_m2_all = 
  p_m2[[1]]+p_m2[[2]]+p_m2[[3]]+p_m2[[4]]+p_m2[[5]]+
  plot_layout(ncol = 2)+
  plot_annotation(
    tag_prefix = 'n = ',
    tag_levels = list(c(100,200,300,500,1000)),
    title = "IPTW estimator mean of survival probabilities under violation of positivity assumption",
  ) &
  theme_classic(8)+
  theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5))
#save as png file
ggsave(p_m2_all,filename="results_plot/figure_survival_probabilities_algorithm2_all.png",width = 9, height = 8) 

##########################
#Empirical standard error#
##########################
make_plot2_se = function(data){
  data%>%
    ggplot(aes(x = Time,y = Empirical_SE,color = factor(tau)))+
    geom_point(aes(shape = factor(tau)))+
    geom_line(aes(linetype = factor(Treatment_regime)))+
    scale_x_continuous(
      breaks = seq(1,5,1),
      labels = seq(1,5,1),
      limits = c(1,5),
      name = 'Time'
    )+
    scale_y_continuous(
      name = "Survival probabilities SE "
    )+
    scale_color_manual(
      values = c("#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
    guides(color=guide_legend(title="Tau"),
           shape = guide_legend(title="Tau"),
           linetype = guide_legend(title="Treatment regime"))
}

# #The list of all plots(p_m2) for different sample size
# p_m2_se = list()
# #create plots for each sample size (n)
# for(i in seq(1,length(n_))){
#   #generate plot for each sample size n
#   plots_ = 
#     results_table2 %>%
#     filter(tau != 'True')%>%
#     #reorder the factor level
#     mutate(tau = fct_relevel(tau,'-0.3','0.3','0.6','1','1.5'))%>%
#     nest(data = -n) %>%
#     filter(n==n_[i])%>%
#     mutate(
#       plots = map(data, make_plot2_se)
#     ) %>%
#     pull(plots) 
#   #entitle the plot
#   p_m2_se[[i]] = 
#     plots_[[1]] +
#     plot_annotation(
#       subtitle = glue("n = {n_[i]}"),
#       title = "IPTW estimator SE of survival probabilities under violation of positivity assumption",
#     ) &
#     theme_classic(8)+
#     theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
#           plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
# }
# p_m2_se
# #save the plots
# for(i in seq(1,length(p_m2_se))){
#   ggsave(p_m2_se[[i]],filename=glue("results_plot/figure_survival_probabilities_algorithm2_{n_[i]}.png"),width = 7, height = 6) 
# }

#Or combine all plots in one plot and save as png file
p_m2_se = list()
for(i in seq(1,length(n_))){
  #generate plot for each sample size n
  plots_ = 
    results_table2 %>%
    filter(tau != 'True')%>%
    #reorder the factor level
    mutate(tau = fct_relevel(tau,'-0.3','0.3','0.6','1','1.5'))%>%
    nest(data = -n) %>%
    filter(n==n_[i])%>%
    mutate(
      plots = map(data, make_plot2_se)
    ) %>%
    pull(plots) 
  #add the plot
  p_m2_se[[i]] = plots_[[1]]
}
#combine the plots in one plot    
p_m2_se_all = 
  p_m2_se[[1]]+p_m2_se[[2]]+p_m2_se[[3]]+p_m2_se[[4]]+p_m2_se[[5]]+
  plot_layout(ncol = 2)+
  plot_annotation(
    tag_prefix = 'n = ',
    tag_levels = list(c(100,200,300,500,1000)),
    title = "IPTW estimator SE of survival probabilities under violation of positivity assumption",
  ) &
  theme_classic(8)+
  theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5))
#save as png file
ggsave(p_m2_se_all,filename="results_plot/figure_survival_probabilities_algorithm2_se_all.png",width = 9, height = 8) 
