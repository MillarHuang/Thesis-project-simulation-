#------------------------------------------------------
#------------------------------------------------------
# sim_results.R
# - Summarises the simulation results to create Tables 1 and 2 and Figures 2 and 3
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#------------------
#time increments: 
#this is a grid of times from 0 to 5 at which the cumulative coefficients and survival probabilities are evaluated
#------------------

t.hor=seq(0,5,0.01)

pos.t1=which(t.hor==1)
pos.t2=which(t.hor==2)
pos.t3=which(t.hor==3)
pos.t4=which(t.hor==4)
pos.t5=which(t.hor==5)

#------------------
#cumulative coefficients at visits k=1,2,3,4,5
#See Table 1 in the paper
#------------------

#----
#true values

est.coef.0.true=coef_0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.A.true=coef_A_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag1.true=coef_Alag1_true[c(pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag2.true=coef_Alag2_true[c(pos.t3,pos.t4,pos.t5)]
est.coef.Alag3.true=coef_Alag3_true[c(pos.t4,pos.t5)]
est.coef.Alag4.true=coef_Alag4_true[pos.t5]

#----
#mean of estimates using MSM-IPTW

est.coef.0.msm=colMeans(coef_0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.A.msm=colMeans(coef_A)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag1.msm=colMeans(coef_Alag1)[c(pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag2.msm=colMeans(coef_Alag2)[c(pos.t3,pos.t4,pos.t5)]
est.coef.Alag3.msm=colMeans(coef_Alag3)[c(pos.t4,pos.t5)]
est.coef.Alag4.msm=colMeans(coef_Alag4)[c(pos.t5)]

#----
#SD of estimates (empirical standard errors) using MSM-IPTW

empsd.coef.0.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_0[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.A.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_A[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag1.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag1[,x])})[c(pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag2.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag2[,x])})[c(pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag3.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag3[,x])})[c(pos.t4,pos.t5)]
empsd.coef.Alag4.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag4[,x])})[c(pos.t5)]

#----
#bias in estimates obtained using MSM-IPTW

bias.coef.0.msm=colMeans(coef_0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-coef_0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.A.msm=colMeans(coef_A)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-coef_A_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.Alag1.msm=colMeans(coef_Alag1)[c(pos.t2,pos.t3,pos.t4,pos.t5)]-coef_Alag1_true[c(pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.Alag2.msm=colMeans(coef_Alag2)[c(pos.t3,pos.t4,pos.t5)]-coef_Alag2_true[c(pos.t3,pos.t4,pos.t5)]
bias.coef.Alag3.msm=colMeans(coef_Alag3)[c(pos.t4,pos.t5)]-coef_Alag3_true[c(pos.t4,pos.t5)]
bias.coef.Alag4.msm=colMeans(coef_Alag4)[c(pos.t5)]-coef_Alag4_true[c(pos.t5)]

#----
#Monte Carlo standard errors for estimates obtained using MSM-IPTW

mcbias.coef.0.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_0[,x]-coef_0_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.A.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_A[,x]-coef_A_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag1.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag1[,x]-coef_Alag1_true[x])/sqrt(n.sim)})[c(pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag2.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag2[,x]-coef_Alag2_true[x])/sqrt(n.sim)})[c(pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag3.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag3[,x]-coef_Alag3_true[x])/sqrt(n.sim)})[c(pos.t4,pos.t5)]
mcbias.coef.Alag4.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag4[,x]-coef_Alag4_true[x])/sqrt(n.sim)})[c(pos.t5)]

#----
#table of results: as in Table 1

#construct the table
table.coef.0=data.frame(Cumulative_coefficient="alpha_0(t)",Time=1:5,
                        True_value=dp3(est.coef.0.true),
                        Mean_estimate__Empirical_SE=table.func(est.coef.0.msm,empsd.coef.0.msm),
                        Bias__MonteCarlo_SE=table.func(bias.coef.0.msm,mcbias.coef.0.msm))

table.coef.A=data.frame(Cumulative_coefficient="alpha_A,0(t)",Time=1:5,
                        True_value=dp3(est.coef.A.true),
                        Mean_estimate__Empirical_SE=table.func(est.coef.A.msm,empsd.coef.A.msm),
                        Bias__MonteCarlo_SE=table.func(bias.coef.A.msm,mcbias.coef.A.msm))

table.coef.Alag1=data.frame(Cumulative_coefficient="alpha_A,1(t)",Time=2:5,
                            True_value=dp3(est.coef.Alag1.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag1.msm,empsd.coef.Alag1.msm)[1:4],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag1.msm,mcbias.coef.Alag1.msm)[1:4])

table.coef.Alag2=data.frame(Cumulative_coefficient="alpha_A,2(t)",Time=3:5,
                            True_value=dp3(est.coef.Alag2.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag2.msm,empsd.coef.Alag2.msm)[1:3],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag2.msm,mcbias.coef.Alag2.msm)[1:3])

table.coef.Alag3=data.frame(Cumulative_coefficient="alpha_A,3(t)",Time=4:5,
                            True_value=dp3(est.coef.Alag3.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag3.msm,empsd.coef.Alag3.msm)[1:2],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag3.msm,mcbias.coef.Alag3.msm)[1:2])

table.coef.Alag4=data.frame(Cumulative_coefficient="alpha_A,4(t)",Time=5,
                            True_value=dp3(est.coef.Alag4.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag4.msm,empsd.coef.Alag4.msm)[1],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag4.msm,mcbias.coef.Alag4.msm)[1])

table1=rbind(table.coef.0,table.coef.A,table.coef.Alag1,table.coef.Alag2,table.coef.Alag3,table.coef.Alag4)
table1$tau = rep(tau,nrow(table1))
table1$n = rep(n,nrow(table1))

#save table as a csv file

write.table(table1,glue("results_data/violation_result_table1_algorithm2_{tau}_{n}.csv"),sep=",",row.names = F)

#------------------
#survival probabilities at visits k=1,2,3,4,5
#See Table 2 in the paper
#------------------

#----
#true values

est.surv1.true=surv1_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.surv0.true=surv0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#mean of estimates using MSM-IPTW

est.surv1.msm=colMeans(surv1)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.surv0.msm=colMeans(surv0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#SD of estimates (empirical standard errors) using MSM-IPTW

empsd.surv1.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv1[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.surv0.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv0[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#bias in estimates obtained using MSM-IPTW

bias.surv1.msm=colMeans(surv1)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-surv1_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.surv0.msm=colMeans(surv0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-surv0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#MC standard errors for estimates obtained using MSM-IPTW

mcbias.surv1.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv1[,x]-surv1_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.surv0.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv0[,x]-surv0_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#latex tables for results: as in Table 2

table.surv.never.treated=data.frame(Treatment_regime="Never treated",Time=1:5,
                                    True_value=dp3(est.surv0.true),
                                    Mean_estimate__Empirical_SE=table.func(est.surv0.msm,empsd.surv0.msm),
                                    Bias__MonteCarlo_SE=table.func(bias.surv0.msm,mcbias.surv0.msm))

table.surv.always.treated=data.frame(Treatment_regime="Always treated",Time=1:5,
                                     True_value=dp3(est.surv1.true),
                                     Mean_estimate__Empirical_SE=table.func(est.surv1.msm,empsd.surv1.msm),
                                     Bias__MonteCarlo_SE=table.func(bias.surv1.msm,mcbias.surv1.msm))

table2=rbind(table.surv.never.treated,table.surv.always.treated)
table2$tau = rep(tau,nrow(table2))
table2$n = rep(n,nrow(table2))
#save table as a csv file

write.table(table2,glue("results_data/violation_result_table2_algorithm2_{tau}_{n}.csv"),sep=",",row.names = F)

