####1.0 LOAD PACKAGES, DATAFRAME, AND PREP####
#packages
library(readxl)
library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MCMCglmm)
library(ggsci)
library(ggpubr)

#data
all_synchrony_data <- read.csv("FILE PATH/all_synchrony_data.csv")
all_synchrony_data <- all_synchrony_data %>%
  mutate_at(vars(species,host,year,culu,sex),as.factor) %>%
  mutate_at(vars(survp,masspm,durationl,day),as.numeric)

####1.1 MODELS####
a <- 1000

#1.1.1 Survival
prior_surv_main<-list(R=list(V=1, fix=1),
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G3=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G4=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G5=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

syncmod_surv_main<-MCMCglmm(survp ~ day + year + day:year, 
                            random = ~ culu +idh(1+day):species + idh(1+day):host + 
                              idh(1+day):species:host + idh(1+day):year:host,
                            family="categorical",
                            prior=prior_surv_main,
                            nitt=55000*200, thin=1*200, burnin=5000*200,
                            verbose=TRUE, pr=TRUE,
                            data=all_synchrony_data)

save(syncmod_surv_main,file="FILE PATH/syncmod_surv_main")

#1.1.2 Mass
#filter data to only species with masspm data, i.e. those which survived to pupation
mass_synchrony_data<-all_synchrony_data[which(all_synchrony_data$masspm > 0),]

#log transform masspm data for analysis
mass_synchrony_data$log_masspm <- log(mass_synchrony_data$masspm)

prior_mass_main<-list(R=list(V=1,nu=0.002),
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G3=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G4=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G5=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G6=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G7=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G8=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

syncmod_mass_main<-MCMCglmm(log_masspm ~ day + sex + year + day:year + day:sex, 
                            random= ~ culu +
                              idh(1+day):species + 
                              idh(1+day):host + 
                              idh(1+day):species:sex + 
                              idh(1+day):host:sex +
                              idh(1+day):species:host + 
                              idh(1+day):year:host +
                              idh(1+day):sex:species:host,
                            family="gaussian",
                            nitt=15000*2000, thin=1*2000, burnin=5000*2000,
                            prior=prior_mass_main, verbose=TRUE,
                            pr=TRUE,
                            data=mass_synchrony_data)

save(syncmod_mass_main,file="FILE PATH/syncmod_mass_main")

#1.1.3 Development time
#filter data to only species with development time data, i.e. those which survived to pupation
dev_synchrony_data<-all_synchrony_data[which(all_synchrony_data$durationl > 0),]

prior_dev_main <-list(R=list(V=1,nu=0.002),
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G2=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G3=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G4=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G5=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G6=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G7=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                             G8=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

syncmod_dev_main<-MCMCglmm(durationl ~ day + sex + year + day:year + day:sex, 
                           random= ~ culu +
                             idh(1+day):species + 
                             idh(1+day):host + 
                             idh(1+day):species:sex + 
                             idh(1+day):host:sex +
                             idh(1+day):species:host + 
                             idh(1+day):year:host +
                             idh(1+day):sex:species:host,
                           family="gaussian",
                           nitt=15000*2000, thin=1*2000, burnin=5000*2000,
                           prior=prior_dev_main, verbose=TRUE,
                           pr=TRUE,
                           data=dev_synchrony_data)

save(syncmod_dev_main,file="FILE PATH/syncmod_dev_main")
