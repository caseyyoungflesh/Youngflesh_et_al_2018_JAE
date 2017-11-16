##############################
#Study R/JAGS code with results
#
#Code for Youngflesh et al. 2017 JAE
#Author: Casey Youngflesh
##############################

#Rethinking ‘normal’: The role of stochasticity in the breeding phenology of a synchronously breeding seabird

#Casey Youngflesh, Stephanie Jenouvrier, Jefferson T. Hinke, Lauren DuBois, Judy St. Leger, Wayne Z. Trivelpiece, Susan G. Trivelpiece, Heather J. Lynch


#Initial set up  - packages

#install packages if they don't exist - then load them
if('pacman' %in% rownames(installed.packages()) == FALSE)
{
  install.packages('pacman')
}

pacman::p_load(ggplot2, 
               reshape2, 
               dplyr, 
               rjags,
               moments,
               MCMCvis)


#Inter-annual variation in CID

#Load data


setwd('Data')

captive_data <- read.csv('Captive_CID.csv', header = TRUE)
wild_data <- read.csv('Wild_CID.csv', header = TRUE)


#Function to calculate intra-annual median and var

srt.fun <- function(IN)
{
  #IN <- SD_lay
  yrs <- unique(IN$YEAR)
  
  OUT <- c()
  for(i in min(yrs):max(yrs))
  {
    #i <- 1993
    temp <- filter(IN, YEAR == i)
    t_md <- median(temp$J_CID)
    t_var <- var(temp$J_CID)
    temp2 <- data.frame(YEAR = i, MEDIAN = t_md, VAR = t_var)
    OUT <- rbind(OUT, temp2)
  }
  return(OUT)
}


#Function to calculate SE of variance
#$\sigma_{s^2} = s^2 * \sqrt{2/(n-1)}$
  
se_var <- function(data) 
{
  OUT <- var(data)*sqrt(2/(length(data)-1))
  return(OUT)
}


#Captive inter-annual variance - $var(y_{i}.)$
captive_md_sd <- srt.fun(captive_data)

var(captive_md_sd$MEDIAN)


#Captive standard error of variance
se_var(captive_md_sd$MEDIAN)


#Wild inter-annual variance - $var(y_{i}.)$
wild_md_sd <- srt.fun(wild_data)

var(wild_md_sd$MEDIAN)


#Wild standard error of variance
se_var(wild_md_sd$MEDIAN)




#Inter-annual variation in CID from literature

#data from Both et al. 2009 and Valtonen et al. 2017

setwd('Data')

lit_data <- read.csv('Lit_data.csv', header = TRUE)

#detrend because we are interested in the interannual variance - trends would conflate this estimate
fit_BT_B_2009 <- lm(lit_data$BT_B_2009 ~ lit_data$YEAR)
res_BT_B_2009 <- residuals(fit_BT_B_2009)
v_BT_B_2009 <- var(res_BT_B_2009)

fit_GT_B_2009 <- lm(lit_data$GT_B_2009 ~ lit_data$YEAR)
res_GT_B_2009 <- residuals(fit_GT_B_2009)
v_GT_B_2009 <- var(res_GT_B_2009)

fit_PF_B_2009 <- lm(lit_data$PF_B_2009 ~ lit_data$YEAR)
res_PF_B_2009 <- residuals(fit_PF_B_2009)
v_PF_B_2009 <- var(res_PF_B_2009)

#no need to detrend as there is no significant change over time
v_CR_V_2017 <- var(lit_data$CR_V_2017, na.rm = TRUE)
v_GT_V_2017 <- var(lit_data$GT_V_2017, na.rm = TRUE)
v_PF_V_2017 <- var(lit_data$PF_V_2017, na.rm = TRUE)







###Plot median CID over time - Captive
plt_t <- data.frame(YEAR = 1992:2015, MD_CID = captive_md_sd$MEDIAN)
plt <- melt(plt_t, id = 'YEAR')


#CAPTIVE PLOT
ggplot(plt, aes(YEAR, value)) +
  geom_line(size = 1.2, col = 'red') +
  theme_bw() +
  ggtitle('Captive penguin breeding phenology') +
  xlab('Year') +
  ylab('CID (days from Sep 30)') +
  coord_cartesian(xlim = c(1990, 2015)) +
  coord_cartesian(ylim = c(25, 45)) +
  scale_x_continuous(breaks = seq(1990, 2015, by = 5)) +
  scale_y_continuous(breaks = c(25, 30, 35, 40, 45)) +
  theme(
    axis.text = element_text(size = 12), #axis label size
    axis.title = element_text(size = 14),
    panel.grid.major = element_line(color = 'gray40'), #lower # is darker
    panel.grid.minor = element_line(color = 'gray95'), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color= 'black', size = 1.5),
    axis.ticks.length= unit(0.15, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1)
  )


###Plot median CID over time - Wild
plt_t2 <- data.frame(YEAR = 1986:2012, MD_CID = wild_md_sd$MEDIAN)
plt2 <- melt(plt_t2, id = 'YEAR')


#WILD PLOT
ggplot(plt2, aes(YEAR, value)) +
  geom_line(size = 1.2, col = 'blue') +
  theme_bw() +
  ggtitle('Wild penguin breeding phenology') +
  xlab('Year') +
  ylab('CID (days from Sep 30)') +
  coord_cartesian(xlim = c(1984, 2012)) +
  coord_cartesian(ylim = c(25, 45)) +
  scale_x_continuous(breaks = seq(1985, 2012, by = 5)) +
  scale_y_continuous(breaks = c(25, 30, 35, 40, 45)) +
  theme(
    axis.text = element_text(size = 12), #axis label size
    axis.title = element_text(size = 14),
    panel.grid.major = element_line(color = 'gray40'), #lower # is darker
    panel.grid.minor = element_line(color = 'gray95'), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color= 'black', size = 1.5),
    axis.ticks.length= unit(0.15, 'cm'), #length of axis tick
    axis.ticks = element_line(size = 1)
  )



#Determine which female is the first to lay in each year (including ties for first)
FEM <- c()
for (i in 1992:2015)
{
  #i <- 1997
  temp <- filter(captive_data, YEAR == i)
  pos <- which(temp$J_CID == min(temp$J_CID))
  t_data <- temp[pos,]
  
  FEM <- rbind(FEM, t_data)
}


#16 different 'leaders' in 24 years

length(unique(FEM$FEMALE_ID))




#Intra-annual variation in CID

#t-test to determine if intra-annual variation differs between captive and wild popualtions
t.test(captive_md_sd$VAR, wild_md_sd$VAR)


#Number of breeders

len_fun <- function(IN)
{
  yrs <- range(IN$YEAR)
  
  LEN <- c()
  for(i in yrs[1]:yrs[2])
  {
    #i <- 1992
    temp <- filter(IN, YEAR == i)
    tl <- dim(temp)[1]
    tb <- c(i, tl)
    LEN <- rbind(LEN, tb)
  }
  
  return(LEN)
}

#range
range(len_fun(captive_data)[,2])


#Hierarchical model - captive population


#JAGS model

#$$y_{ij} = \mu + \alpha_{i} + \beta_{j} + \gamma*AGE_{ij} + \epsilon_{ij}$$
#$$\alpha_{i} \sim N(0, \tau_{year})$$
#$$\beta_{j} \sim N(0, \tau_{individual})$$
#$$\epsilon_{ij} \sim N(0, \tau_{model})$$

setwd('Data')

AGE_mat <- read.csv('AGE_mat.csv', header= TRUE)
CID_mat <- read.csv('CID_mat.csv', header= TRUE)

DATA <- list(
  y = CID_mat,
  yr = 1:NCOL(CID_mat),
  ind = 1:NROW(CID_mat),
  age = AGE_mat,
  N = NCOL(CID_mat), #columns are year in matrix
  M = NROW(CID_mat)) #rows are individuals


#----------------#
#model

#alpha = YEAR - random
#beta = INDIVIDUAL - random
#gamma = AGE - fixed

setwd('..')


{
  sink("captive.jags")
  
  cat("
      model {
      
      for(t in 1:N)
      {
      for(i in 1:M)
      {
      y[i,t] ~ dnorm(mu.g[i,t], tau)
      mu.g[i,t] <- mu + alpha[yr[t]] + beta[ind[i]] + gamma*age[i,t]
      }
      }
      
      #priors
      
      #year
      for(t in 1:N)
      {
      alpha[t] ~ dnorm(0, tau.year)
      }
      
      #individual
      for(i in 1:M)
      {
      beta[i] ~ dnorm(0, tau.ind)
      }
      
      #mu, gamma, and tau
      mu ~ dnorm(0, 0.001)
      gamma ~ dnorm(0, 0.001)
      tau ~ dgamma(0.01, 0.01)    
      
      #hyperparameters
      tau.year ~ dgamma(0.01, 0.01)
      tau.ind ~ dgamma(0.01, 0.01)
      rate.tau ~ dgamma(0.01, 0.01)
      shape.tau ~ dgamma(0.01, 0.01)
      
      }",fill = TRUE)

  sink()
}


#Run model

#----------------------#
#Starting values


Inits <- function() {list(alpha = rep(rnorm(1), 
                                      ncol(CID_mat)),
                          beta = rep(rnorm(1), 
                                     nrow(CID_mat)), 
                          tau = rgamma(1, 1),
                          mu = rnorm(1), 
                          gamma = rnorm(1), 
                          tau.year = rgamma(1, 1), 
                          tau.ind = rgamma(1, 1),
                          rate.tau = rgamma(1, 1), 
                          shape.tau = rgamma(1, 1))}


#----------------------#
#Parameters to track

Pars <- c("alpha", 
          "beta", 
          "gamma")


# Inputs for MCMC ---------------------------------------------------------

n_adapt <- 5000  # number for initial adapt
n_burn <- 40000 # number burnin
n_draw <- 10000  # number of final draws to make
n_thin <- 2    # thinning rate
n_chain <- 3  # number of chains


#----------------------#
#Run model

jm = jags.model(data = DATA, 
                file = "captive.jags", 
                inits = Inits, 
                n.chains = 3, 
                n.adapt = n_adapt)

update(jm, n.iter = n_burn)

out = coda.samples(jm, 
                   n.iter = n_draw, 
                   variable.names = Pars, 
                   thin = n_thin)


###Year effect (alpha)

MCMCplot(out, 
         params = 'alpha', 
         rank = FALSE, 
         labels = NULL,
         horiz = FALSE, 
         ref_ovl = FALSE, 
         ylim = c(-15, 15))

MCMCsummary(out, 
            params = 'alpha')


###Individual effect (beta)

MCMCplot(out, 
         params = 'beta', 
         rank = TRUE, 
         labels = NULL,
         horiz = FALSE, 
         thick_sz = 2, 
         thin_sz = 1, 
         med_sz = .6,
         ref_ovl = FALSE, 
         ylim = c(-15, 15))

MCMCsummary(out, 
            params = 'beta')

###Age effect (gamma)

MCMCplot(out, 
         params = 'gamma', 
         rank = TRUE, 
         labels = NULL,
         horiz = FALSE, 
         thick_sz = 10, 
         thin_sz = 3, 
         med_sz = 3,
         ref_ovl = FALSE, 
         ylim = c(-15, 15))

MCMCsummary(out, 
            params = 'gamma')


#Skewness test and plots

###Captive

#scale and aggregate data

sc_agg_fun <- function(IN)
{
  yrs <- range(IN$YEAR)
  
  OUT <- c()
  for(i in yrs[1]:yrs[2])
  {
    temp <- filter(IN, YEAR == i)
    s_CID <- scale(temp$J_CID, scale = TRUE)
    
    t.out <- cbind(temp, s_CID)
    OUT <- rbind(OUT, t.out)
  }
  return(OUT)
}


cap_sk <- sc_agg_fun(captive_data)


#Skew determined using D'Agostino test
agostino.test(cap_sk$s_CID)


#Standard error skew
#standard error of skewness function
ses <- function(n) 
{
  sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)))
}

length_captive <- dim(cap_sk)[1]

ses(length_captive)


###Wild
wild_sk <- sc_agg_fun(wild_data)


#Skew
agostino.test(wild_sk$s_CID)

#Standard error skew
length_wild <- dim(wild_sk)[1]

ses(length_wild)



###Simulate normal breeding distribution given true CID mean and sd - Plot 

#Captive
m_cap <- mean(cap_sk$s_CID, na.rm = TRUE)
sd_cap <- sd(cap_sk$s_CID, na.rm = TRUE)
cap_rd <- rnorm(100000, mean = m_cap, sd = sd_cap)

hist(cap_sk$s_CID, prob = TRUE, 
main = 'Breeding distribution - Captive',
xlab = 'CID', ylab= 'Density', col = 'grey90', 
xlim = c(-2.5, 3),
breaks = 15)
lines(density(cap_rd), col = 'blue', lwd = 5)
lines(density(cap_sk$s_CID, na.rm = TRUE), col = 'red', lwd = 5)


#Wild
m_wild <- mean(wild_sk$s_CID, na.rm=TRUE)
sd_wild <- sd(wild_sk$s_CID, na.rm=TRUE)
wild_rd <- rnorm(100000, mean= m_wild, sd= sd_wild)

hist(wild_sk$s_CID, prob=TRUE, 
main='Breeding distribution - Wild',
xlab= 'CID', ylab= 'Density', col = 'grey90', 
xlim = c(-2.5, 3), 
breaks = 25)
lines(density(wild_rd), col='blue', lwd=5)
lines(density(wild_sk$s_CID, na.rm=TRUE), col='red', lwd=5)









###Skew estimates from literature


setwd('Data')

lit_skew <- read.csv('Lit_skew.csv', header = TRUE)

mn_skew_lit  <- mean(lit_skew$Skew)

hist(lit_skew$Skew, 
     main = 'Red = Captive; Blue = Wild',
     xlab = 'Estimated Skew')
abline(v = skew_captive, col = 'red', lwd = 5)
abline(v = skew_wild, col = 'blue', lwd = 5)


per_captive <- length(which(lit_skew$Skew > skew_captive))/NROW(lit_skew)
per_wild <- length(which(lit_skew$Skew > skew_wild))/NROW(lit_skew)






###Simulation of right skew breeding distribution


#simulate a normal distirbution (breeding in the absence of conspecifics)
set.seed(1)
time <- round(rnorm(1000, 50, 10))
time <- sort(time)
#no significant skew
agostino.test(time)

#simulate 'contagion effect' of breeding
for (i in min(time):max(time))
{
  time[time == i] <- time[time == i] - sum(as.numeric((time > (i - 3)) & (time < i)))/4
}

hist(time, prob = TRUE)
lines(density(time), col = 'blue', lwd = 5)

agostino.test(time)

