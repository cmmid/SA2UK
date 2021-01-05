# - - 
# Plot introductions vs growth

# Helper functions

library(tidyverse)

c.text <- function(x,sigF=3){
  if(is.numeric(x)){bp1=signif(x,sigF)}else{bp1=x}
  paste(bp1[1]," (95% CI:",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}


# Baseline parameters
R <- 1.5 # Local reproduction number

gt <- 6.8 # Generation time
lambda_UK <- R/gt # lambda for transmission model

detected_cases <- 2 # Initial reported cases
detection_date <- as.Date("2020-12-12")
  
prop_inf_rep_ONS <- 0.6
prop_inf_rep_REACT <- 0.3
prop_inf_rep_mean <- 0.45

prop_sequenced <- 0.057

# Estimate total cases

# Function to calculate CI for total cases
find_CI <- function(detected_cases,pp){
  range_nn <- seq(detected_cases,1e4*detected_cases,1)
  lik <- dbinom(detected_cases,range_nn,pp,log=T)
  
  lik_surface <- lik - max(lik) # Find max likelihood
  
  CI_range <- range_nn[lik_surface>-1.92] # Find possible values
  
  c(round(detected_cases/pp),min(CI_range),max(CI_range))
  
}

# Epidemic size

final_size <- function(ss,RR){abs(RR*(ss-1)-log(ss))}

attack_rate_base <- 100*(1-optimize(final_size, interval = c(0,1),RR=1.2)$minimum)
attack_rate_50 <- 100*(1-optimize(final_size, interval = c(0,1),RR=1.2*1.5)$minimum)


# 
ss_range <- seq(0,1,0.01)

final_size(ss_range,RR=1.5)

plot(ss_range,final_size(ss_range,RR=2))

#

total_US <- find_CI(1,pp= 0.5 * 0.0025)

# Total cases

total_SA_variant_inf_ONS <- find_CI(detected_cases,pp=prop_inf_rep_ONS* prop_sequenced)
total_SA_variant_inf_REACT <- find_CI(detected_cases,pp=prop_inf_rep_REACT* prop_sequenced)
total_SA_variant_inf_mean <- find_CI(detected_cases,pp=prop_inf_rep_mean* prop_sequenced)

total_SA_variant_inf_mean %>% c.text()


# Need solve: I(t) = a/lambda * ( exp(lambda*t) - 1 ) = total_SA_variant_infections

# Estimate time since first introduction
time_intro <- function(lambda,a,total_n){log(total_n*(lambda/a) + 1)/lambda}

aa1 <- 2
aa2 <- 4

yy_2 <- time_intro(lambda_UK,a=aa1,total_SA_variant_inf_mean) # Estimate import time
yy_4 <- time_intro(lambda_UK,a=aa2,total_SA_variant_inf_mean) # Estimate import time

yy_2 %>% c.text(.,sigF=2)
yy_4 %>% c.text(.,sigF=2)

c((detection_date - yy_2)[3], (detection_date - yy_4)[2] ) # Range of dates of introduction


# Estimate total imported cases up to most recent detection date

(aa1*yy_2) %>% c.text(.,sigF=2)
(aa2*yy_4) %>% c.text(.,sigF=2)

# par(mfrow=c(1,2),mar=c(3.5,3.5,1,1),mgp=c(2,0.6,0),las=1)
# 
# plot(aa,yy,type="l",xlab="daily imported infections",ylab="time since first import",xaxs="i",yaxs="i",ylim=c(0,20))
# lines(aa,yy1,col="blue")
# 
# plot(aa,yy*aa,type="l",xlab="daily imported infections",ylab="estimated total imports",xaxs="i",yaxs="i",ylim=c(0,100))
# lines(aa,yy1*aa,col="blue")




