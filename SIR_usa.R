rm(list = ls())

setwd('/Users/chrissymo/Documents/MSIS/research/with Xialu/Case and death counts fetch(data for Rt)')

library(dplyr)
library(ggplot2)
library(epidemia)
library(rstanarm)
#library(readxl)
library(car)
#import data
#jhu_complete <- read_excel("Documents/MSIS/research/with Xialu/Case and death counts fetch/jhu_complete_v4.xlsx")
#data_test <- jhu_complete_v4
#data_test <- jhu_mainstates_1to11_v3
completed_data <- jhu_mainstates_1to11_v3
#main_states <- jhu_mainstates_1to11_v3
data_test <- completed_data
#data_test <- jhu_complete_1to11
#head(data_test)
#data_test <- filter(data_test,date > date[which(cumsum(deaths) > 10)[1]-30])
data_test <- data_test%>%group_by(state) %>% mutate(csum = cumsum(deaths))
data_test <- filter(data_test,date > date[which(csum > 10)[1]-30])
summarise(data_test, start = min(date), end = max(date))
  
  
##model components
#transmission model
rt_test <- epirt(
    formula = R(state, date) ~ 1 + lockdown + reopen,
    prior = shifted_gamma(shape=1/6, scale = 1, shift = log(1.05)/6),
    prior_covariance = decov(shape = c(2, rep(0.5, 5)),scale=0.25),
    link = scaled_logit(6.5)
  )
  #infections
inf_test <- epiinf(gen = EuropeCovid2$si, seed_days = 6)
  
  #observations
deaths_test <- epiobs(formula = deaths ~ 1, i2o = EuropeCovid2$inf2death,
                   prior_intercept = normal(0,0.2), link = scaled_logit(0.02))

#########################################
##work on one state
data <- filter(data_test, state == 'wa')
#data <- data_test
args_test <- list(rt = rt_test, inf = inf_test, obs = deaths_test, data = data, seed = 12345,
                  refresh = 0)
##=Approximating the Posterior
args_test$algorithm <- "fullrank"; args_test$iter <- 5e4; args_test$tol_rel_obj <- 1e-3
fm_test<- do.call(epim, args_test)
#plot rt
plot_rt(fm_test, step = T, levels = c(50,95)) + theme_bw()

rt_value <- posterior_rt(fm_test)
true_rt = apply(rt_value$draws, 2, function(x) quantile(x, 0.5))
unique(true_rt)

#plot death
plot_obs(fm_test, type = "deaths", levels = c(50, 95))

#plot infections
plot_infections(fm_test)
plot_infectious(fm_test)


#########################################
###work on all states

states = unique(data_test$state)

rt_values_before <- c()
rt_values_after <- c()
rt_values_reopen <- c()
for (i in 1:length(states)) {
    ##assign data a state or whole states
  data <- filter(data_test, state == states[i])
    
    ##model fitting
    #prior check
  args_test <- list(rt = rt_test, inf = inf_test, obs = deaths_test, data = data, seed = 12345,
                 refresh = 0)
    #pr_args_test <- c(args_test, list(algorithm = "sampling", iter = 1e3, prior_PD = TRUE))
    #fm_prior_test <- do.call(epim, pr_args_test)
    #plot rt
    #plot_rt(fm_prior_test) + theme_bw()
    
    ##=Approximating the Posterior
  args_test$algorithm <- "fullrank"; args_test$iter <- 5e4; args_test$tol_rel_obj <- 1e-3
  fm_test<- do.call(epim, args_test)
    #plot rt
    #plot_rt(fm_test, step = T, levels = c(50,95)) + theme_bw()
    
  rt_value <- posterior_rt(fm_test)
  true_rt = apply(rt_value$draws, 2, function(x) quantile(x, 0.5))
  unique_true_rt = unique(true_rt)

  rt_values_before[i] <- unique_true_rt[1]
  rt_values_after[i] <- unique_true_rt[2]
  rt_values_reopen[i] <- unique_true_rt[3]
}

#df_rt <- data.frame(states,rt_values_before,rt_values_after)
df_rt <- data.frame(states,rt_values_before,rt_values_after,rt_values_reopen)
#save results
write.csv(df_rt,'df_rt_add_reopen_1105.csv')



library(MASS)
data_lm <- input_data_for_regression_10_19
data_lm <- data_lm[-c(7,49,91),] ##delete DC
#data_lm <- data_lm[-c(29,21),]##delete outliners

#regular regression
lm_re_data <- subset(data_lm,select = -c(id,state))
lm_re_data <- subset(data_lm,select = -c(id,state,Households_Median_income_dollars,Percentage_bachelor_degree_or_higher))
#lm_re_data <- subset(data_lm,select = c(Rt,lock_down,GDP_per_capita, Density))
lm_re_model <- lm(formula = Rt~.,data = lm_re_data)
lm_re_model <- lm(formula = log(Rt)~.,data = lm_re_data)
lm_re_model <- lm(formula = log(Rt)~lock_down+reopen+Density+
                    log(GDP_per_capita)+log(Households_Median_income_dollars)+
                    Top_enplaned_or_not+Republican,data = lm_re_data)


#add interaction between lock down and density
lm_re_data$lock_density <- lm_re_data$lock_down * lm_re_data$Density
lm_re_model <- lm(log(Rt)~lock_down + reopen + Density + GDP_per_capita + lock_density + lock_gdp, data = lm_re_data)
lm_re_model <- lm(formula = log(Rt)~lock_down+reopen+Density+log(GDP_per_capita)
                  +lock_down*Density+reopen*Density,data = lm_re_data)

#summary model
summary(lm_re_model)
vif(lm_re_model)
qqPlot(lm_re_model$residuals) 

#interaction regression
#data_lm = interaction_factors
#lm_interaction_model <- lm(formula = Rt~.,data = data_lm)
#summary(lm_interaction_model)

#using boxcox to find out lambda
bc <- boxcox(lm_re_model)
lambda <- bc$x[which.max(bc$y)]
new_model <- lm(((Rt^lambda-1)/lambda) ~ ., data = lm_re_data)
new_model <- lm(((Rt^lambda-1)/lambda) ~ lock_down + Density + GDP_per_capita + lock_density + lock_gdp, data = lm_re_data)
summary(new_model)
qqPlot(new_model$residuals) 



##plot cases & deaths
path <-'/Users/chrissymo/Documents/MSIS/research/with Xialu/Case and death counts fetch/plot/'

for (i in 1:length(states)) {
  data <- filter(data_test, state == states[i])
  plot_name <- paste(toString(states[i]),"deaths")
  jpeg(paste(path,plot_name,'.jpg'))
  plot(data$date,data$deaths,type='l',main = plot_name, xlab = 'date', ylab ='deaths')
  dev.off()
}





############## 52 states Rt plot ##############################
path <-'/Users/chrissymo/Documents/MSIS/research/with Xialu/Case and death counts fetch(data for Rt)/52states_plots/'
states = unique(data_test$state)
#length(states)
for (i in length(states)) {
  data <- filter(data_test, state == states[i])
  plot_name <- paste(toString(states[i]),"Rt")
  jpeg(paste(path,plot_name,'.jpg'))
  #epim component
  args_test <- list(rt = rt_test, inf = inf_test, obs = deaths_test, data = data, seed = 12345,
                    refresh = 0)
  ##Approximating the Posterior
  args_test$algorithm <- "fullrank"; args_test$iter <- 5e4; args_test$tol_rel_obj <- 1e-3
  fm_test<- do.call(epim, args_test)
  plot_rt(fm_test, step = T, levels = c(50,95)) + theme_bw()
  dev.off()
}



data <- filter(data_test, state == states[1])
plot_name <- paste(toString(states[1]),"Rt")
jpeg(paste(path,plot_name,'.jpg'))
#epim component
args_test <- list(rt = rt_test, inf = inf_test, obs = deaths_test, data = data, seed = 12345,
                  refresh = 0)
##Approximating the Posterior
args_test$algorithm <- "fullrank"; args_test$iter <- 5e4; args_test$tol_rel_obj <- 1e-3
fm_test<- do.call(epim, args_test)
#plot rt
plot_rt(fm_test, step = T, levels = c(50,95)) + theme_bw()
dev.off()