# sim_results_violation_algorithm1.R
# - obtain the results of IPTW estimator mean, standard deviation, bias and MC standard deviation of bias, along different CD4 threshold and different sample size (n)
# - generate the corresponding plots for mean & standard deviation of IPTW estimator under different sample size, combine all plots together as one plot, then save as png file

#Sample size (for each simulation)
n_ = c(100,200,300,500,1000)

#Threshold of CD4 : when subject's CD4 <= tau, the subject is forced to start the treatment 
tau_ = c(-100,100,200,300,400,500)  #tau=-100 means no threshold(the estimate under positivity assumption)
 
#IPTW estimator results on simulated data
IPTW_df = data.frame()
#True values of parameters
True_parameter = c(0.05,-1.5,0.1)
#Simulating data under certain violation threshold(tau) with n sample size, and attain its IPTW estimator
for(n in n_){
  for (tau in tau_){
  #simulate and estimate the IPTW estimator in simulation study with B = 100 repetitions
  source('simulation_analysis_algorithm1.R') 
  df = data.frame(
    "sample size" = n,
    "tau" = tau,
    'gamma1' = round(mean(gamma1),3),
    'gamma1_std' = round(sd(gamma1),3),
    'gamma1_bias' = round(mean(gamma1)-True_parameter[1],3),
    'gamma1_MCstd' = round(sd(gamma1-True_parameter[1])/sqrt(B),3),
    'gamma2' = round(mean(gamma2),3),
    'gamma2_std' = round(sd(gamma2),3),
    'gamma2_bias' = round(mean(gamma2)-True_parameter[2],3),
    'gamma2_MCstd' = round(sd(gamma2-True_parameter[2])/sqrt(B),3),
    'gamma3' = round(mean(gamma3),3),
    'gamma3_std' = round(sd(gamma3),3),
    'gamma3_bias' = round(mean(gamma3)-True_parameter[3],3),
    'gamma3_MCstd' = round(sd(gamma3-True_parameter[3])/sqrt(B),3)
  )
  IPTW_df = rbind(IPTW_df, df)
  }
  print(n)
  # save results of certain sample size n as a csv file
  # write.table(IPTW_df_,glue("results_data/violation_result_algorithm1_{n}.csv"),sep=",",row.names = F)
}

#save results of all sample size as a csv file
write.table(IPTW_df,"results_data/violation_result_algorithm1_all.csv",sep=",",row.names = F)

#import results
results_wide = read_csv("results_data/violation_result_algorithm1_all.csv")
#transform the results to long-format
gamma.long = results_wide %>% 
  select(sample.size,tau,gamma1,gamma2,gamma3)%>%
  pivot_longer(cols = -c(sample.size,tau), names_to='gamma', values_to = 'mean')
gamma_std.long = results_wide %>% 
  select(sample.size,tau,gamma1_std,gamma2_std,gamma3_std)%>%
  pivot_longer(cols = -c(sample.size,tau), names_to='gamma', values_to = 'std')

#Visualization
#function to generate plot for IPTW estimator mean
make_plot_mean = function(data){
  data%>%
    ggplot(aes(x = tau,y = mean,color = factor(sample.size)))+
    geom_point(aes(shape = factor(sample.size)))+
    geom_line(aes(linetype = factor(sample.size)))+
    scale_x_continuous(
      expand = c(0, 50),
      breaks = seq(-100,500,100),
      labels= c('True value',seq(0,500,100)),
      name = 'tau'
    )+
    scale_y_continuous(
      name = glue("Monte Carlo mean"),
    )+
    scale_color_manual(
      values = c("#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
    guides(color=guide_legend(title="Sample size"),
           shape = guide_legend(title="Sample size"),
           linetype = guide_legend(title="Sample size"))
}
#function to generate IPTW empirical standard deviation plot
make_plot_std = function(data){
  gamma = str_split(gamma,"_",simplify = TRUE)[1]
  data%>%
    ggplot(aes(x = tau,y = std,color = factor(sample.size)))+
    geom_point(aes(shape = factor(sample.size)))+
    geom_line(aes(linetype = factor(sample.size)))+
    scale_x_continuous(
      expand = c(0, 50),
      breaks = seq(-100,500,100),
      labels= c('True value',seq(0,500,100)),
      name = 'tau'
    )+
    scale_y_continuous(
      name = glue("Empirical standard deviation"),
    )+
    scale_color_manual(
      values = c("#F8766D","#B79F00" ,"#00BA38", "#00BFC4", "#619CFF" )
    )+
    guides(color=guide_legend(title="Sample size"),
           shape = guide_legend(title="Sample size"),
           linetype = guide_legend(title="Sample size"))
}

gamma_ = c('gamma1','gamma2','gamma3')
#The list of all plots (includes all gamma)
p_m = list()
#create plots of mean and standard deviaiton of IPTW estimator for each gamma
for(i in seq(1,length(gamma_))){
  #generate plot for each cumulative coefficient under a certain sample size n
  plots_mean = 
    gamma.long %>%
    nest(data = -gamma) %>%
    filter(gamma==gamma_[i])%>%
    mutate(
      plots = map(data, make_plot_mean)
    ) %>%
    pull(plots) 
  plots_std = 
    gamma_std.long %>%
    nest(data = -gamma) %>%
    filter(gamma==glue('gamma{i}_std'))%>%
    mutate(
      plots = map(data, make_plot_std)
    ) %>%
    pull(plots) 
  #combine the plots of different cumulative coefficient together for each certain sample size n 
  p_m[[i]] = plots_mean[[1]] | plots_std[[1]]
}
#combine all plots together as one plot
p_all = p_m[[1]]/p_m[[2]]/p_m[[3]]+ 
  plot_annotation(
    tag_levels = list(c('gamma1', '','gamma2','','gamma3','')),
    title = "IPTW estimator mean/standard deviation under violation of positivity assumption",
  ) &
  theme_classic(8)+
  theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))

#save the plot as png file
ggsave(p_all,filename="results_plot/figure_mean_std_algorithm1_all.png",width = 7, height = 6) 



  
# #visualization
# #gamma1
# p1 =gamma.long%>%
#   filter(gamma=='gamma1')%>%
#   ggplot(aes(x = tau, y = mean)) +
#   geom_line(aes(color=gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Monte Carlo mean",
#     )+
#   scale_color_manual(
#     values = '#F85B52' 
#   )
# 
# p1_v =gamma_std.long%>%
#   filter(gamma=='gamma1_std')%>%
#   ggplot(aes(x = tau, y = std)) +
#   geom_line(aes(color=gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Empirical standard deviation",
#   )+
#   scale_color_manual(
#     values = '#F85B52',
#     labels=c('gamma1')
#   )
# #gamma2
# p2 =gamma.long%>%
#   filter(gamma=='gamma2')%>%
#   ggplot(aes(x = tau, y = mean)) +
#   geom_line(aes(color = gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Monte Carlo mean",
#   )+
#   scale_color_manual(
#     values = "#7194BA"
#   )
# 
# p2_v =gamma_std.long%>%
#   filter(gamma=='gamma2_std')%>%
#   ggplot(aes(x = tau, y = std)) +
#   geom_line(aes(color=gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Empirical standard deviation",
#   )+
#   scale_color_manual(
#     values = "#7194BA",
#     labels=c('gamma2')
#   )
# #gamma3
# p3 =gamma.long%>%
#   filter(gamma=='gamma3')%>%
#   ggplot(aes(x = tau, y = mean)) +
#   geom_line(aes(color = gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Monte Carlo mean",
#   )+
#   scale_color_manual(
#     values = "#5C9E76"
#   )
# 
# p3_v =gamma_std.long%>%
#   filter(gamma=='gamma3_std')%>%
#   ggplot(aes(x = tau, y = std)) +
#   geom_line(aes(color=gamma)) +
#   scale_x_continuous(
#     expand = c(0, 50),
#     breaks = seq(-100,500,100),
#     labels= c('True value',seq(0,500,100)),
#     name = 'tau'
#   )+
#   scale_y_continuous(
#     name = "Empirical standard deviation",
#   )+
#   scale_color_manual(
#     values = "#5C9E76",
#     labels=c('gamma3')
#   )
# #Patch plotting
# #Monte Carlo mean
# p_m = p1 / p2 / p3  +
#   plot_annotation(
#     tag_levels = 'A',
#     subtitle = glue("n = {n}"),
#     title = "IPTW estimator mean under violation of positivity assumption",
#   ) &
#    theme_classic(8)+
#    theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
#          plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
# #Empirical standard deviation
# p_std = p1_v / p2_v / p3_v  +
#   plot_annotation(
#     tag_levels = 'A',
#     subtitle = glue("n = {n}"),
#     title = "IPTW estimator standard deviation under violation of positivity assumption",
#   ) &
#   theme_classic(8)+
#   theme(plot.title = element_text(size= 10,face = 'bold', hjust = 0.5),
#         plot.subtitle = element_text(size= 10,face = 'bold', hjust = 0.05))
# 
# #save plot as png file
# ggsave(p_m,filename=glue("results_plot/figure_mean_algorithm1_{n}.png"),width = 7, height = 6) 
# ggsave(p_std,filename=glue("results_plot/figure_std_algorithm1_{n}.png"),width = 7, height = 6) 


