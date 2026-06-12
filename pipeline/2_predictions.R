####2.0 MODEL PREDICTIONS####
#predictions are made across various levels, including (1) the overall average, 
#(2) average across caterpillar species, (3) average across host-plant species,
#(4) trend for each host/caterpillar species pairwise combination

#load required data: all_synchrony_data.csv

a<-1000

####2.1 SURVIVAL####
#generate blank dataframe for predictions
all_synchrony_data_2020<-all_synchrony_data[which(all_synchrony_data$year=="2020"),]
all_synchrony_data_2021<-all_synchrony_data[which(all_synchrony_data$year=="2021"),]
blank_surv_20_b <- as.data.frame(cbind(rep(c("brumata"),80),rep(c("Betula","Crataegus","Malus","Quercus"),each=20),rep((seq(0, 19, 1)),4),rep(c(2020),80),rep(c(1,2,3,4),20),rep(c(0,1,0,1),20)))
names(blank_surv_20_b)<-c("species","host","day","year","culu","survp")
blank_surv_20_d <- as.data.frame(cbind(rep(c("defoliaria"),80),rep(c("Betula","Crataegus","Malus","Quercus"),each=20),rep((seq(0, 19, 1)),4),rep(c(2020),80),rep(c(1,2,3,4),20),rep(c(0,1,0,1),20)))
names(blank_surv_20_d)<-c("species","host","day","year","culu","survp")
blank_surv_20_a <- as.data.frame(cbind(rep(c("antiqua"),531),rep(levels(all_synchrony_data_2020$host),each=59),rep((seq(0, 58, 1)),9),rep(c(2020),531),rep(c(1,2,3),177),rep(c(0,1,0),177)))
names(blank_surv_20_a)<-c("species","host","day","year","culu","survp")
blank_surv_21 <- as.data.frame(cbind(rep(levels(all_synchrony_data_2021$species),each=528),
                                     rep(c("Acer","Alnus","Betula","Crataegus","Malus","Quercus","Scapraea","Salba"),each=66,times=6),
                                     rep((seq(0, 65, 1)),48),
                                     rep(c(2021),3168),
                                     rep(c(1,2,3),1056),
                                     rep(c(0,1,0),1056)))
names(blank_surv_21)<-c("species","host","day","year","culu","survp")
blank_surv<-rbind(blank_surv_20_b, blank_surv_20_d, blank_surv_20_a, blank_surv_21)
blank_surv <- blank_surv %>% mutate_at(c("species","host","year","culu"), as.factor)
blank_surv <- blank_surv %>% mutate_at(c("survp","day"), as.numeric)

####2.1.1 Main trend####
prior_surv_main<-list(R=list(V=1, fix=1),
                      G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G2=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                             G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))
surv_pred_main<-predict(syncmod_surv_main,newdat=blank_surv,marginal=~culu + 
                          idh(1+day):species + 
                          idh(1+day):host + 
                          idh(1+day):species:host + 
                          idh(1+day):year:host,
                        type="response",interval="confidence",prior=prior_surv_main)
surv_predictions_main<-cbind(blank_surv, surv_pred_main)
#cut down to antiqua, because longest 2020 day run, and keep one for each year
surv_predictions_main<-surv_predictions_main[which(surv_predictions_main$species=="antiqua" & surv_predictions_main$host=="Betula"),]
#drop uneccessary rows
surv_predictions_main$culu<-NULL
surv_predictions_main$survp<-NULL
surv_predictions_main$species<-NULL
surv_predictions_main$host<-NULL
write.csv(surv_predictions_main,"FILE PATH/surv_predictions_main.csv", row.names = FALSE)


####2.1.2 Caterpillar species averages####
prior_surv_species<-list(R=list(V=1, fix=1),
                         G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                                G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                                G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                                G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))
surv_pred_species<-predict(syncmod_surv_main,newdat=blank_surv,marginal=~culu + 
                             idh(1+day):host + 
                             idh(1+day):species:host + 
                             idh(1+day):year:host,
                           type="response",interval="confidence", prior=prior_surv_species)
surv_predictions_species<-cbind(blank_surv, surv_pred_species)
#cut down to one set per species
surv_predictions_species<-surv_predictions_species[which(surv_predictions_species$host=="Betula"),]
surv_predictions_species$culu<-NULL
surv_predictions_species$survp<-NULL
surv_predictions_species$host<-NULL
write.csv(surv_predictions_species,"FILE PATH/surv_predictions_species.csv", row.names = FALSE)

####2.1.3 Host species averages####
prior_surv_hosts<-list(R=list(V=1, fix=1),
                       G=list(G1=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                              G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                              G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a),
                              G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*a)))
surv_pred_hosts<-predict(syncmod_surv_main,newdat=blank_surv,marginal=~culu + 
                           idh(1+day):species + 
                           idh(1+day):species:host + 
                           idh(1+day):year:host,
                         type="response",interval="confidence", prior=prior_surv_hosts)
surv_predictions_hosts<-cbind(blank_surv, surv_pred_hosts)
#cut down to one set per host
surv_predictions_hosts<-surv_predictions_hosts[which(surv_predictions_hosts$species=="antiqua"),]
surv_predictions_hosts$culu<-NULL
surv_predictions_hosts$survp<-NULL
surv_predictions_hosts$species<-NULL
write.csv(surv_predictions_hosts,"FILE PATH/surv_predictions_hosts.csv", row.names = FALSE)

####2.1.4 All species combinations####
surv_pred_all<-predict(syncmod_surv_main,newdat=blank_surv,marginal=~culu,type="response",interval="confidence")
surv_predictions_all<-cbind(blank_surv, surv_pred)
surv_predictions_all$culu<-NULL
surv_predictions_all$survp<-NULL
write.csv(surv_predictions_all,"FILE PATH/surv_predictions_all.csv", row.names = FALSE)


####2.2 MASS####
#generate blank dataframe for predictions
blank_mass_20_b <- as.data.frame(cbind(rep(c("brumata"),80),rep(c("Betula","Crataegus","Malus","Quercus"),each=20),rep((seq(0, 19, 1)),4),rep(c(2020),80),rep(c(1,2,3,4),20),rep(c(1,2,3,4),20),rep(c("F"),80)))
names(blank_mass_20_b)<-c("species","host","day","year","culu","log_masspm","sex")
blank_mass_20_d <- as.data.frame(cbind(rep(c("defoliaria"),80),rep(c("Betula","Crataegus","Malus","Quercus"),each=20),rep((seq(0, 19, 1)),4),rep(c(2020),80),rep(c(1,2,3,4),20),rep(c(1,2,3,4),20),rep(c("F"),80)))
names(blank_mass_20_d)<-c("species","host","day","year","culu","log_masspm","sex")
blank_mass_20_a <- as.data.frame(cbind(rep(c("antiqua"),531),rep(levels(all_synchrony_data_2020$host),each=59),rep((seq(0, 58, 1)),9),rep(c(2020),531),rep(c(1,2,3),177),rep(c(1,2,3),177),rep(c("F"),531)))
names(blank_mass_20_a)<-c("species","host","day","year","culu","log_masspm","sex")
blank_mass_21 <- as.data.frame(cbind(rep(levels(all_synchrony_data_2021$species),each=528),
                                     rep(c("Acer","Alnus","Betula","Crataegus","Malus","Quercus","Scapraea","Salba"),each=66,times=6),
                                     rep((seq(0, 65, 1)),48),
                                     rep(c(2021),3168),
                                     rep(c(1,2,3),1056),
                                     rep(c(1,2,3),1056),
                                     rep(c("F"),3168)))
names(blank_mass_21)<-c("species","host","day","year","culu","log_masspm","sex")
blank_mass_F<-rbind(blank_mass_20_b, blank_mass_20_d, blank_mass_20_a, blank_mass_21)

blank_mass_M<-blank_mass_F
blank_mass_M$sex<-"M"

blank_mass<-rbind(blank_mass_F, blank_mass_M)

blank_mass <- blank_mass %>% mutate_at(c("species","host","year","culu","sex"), as.factor)
blank_mass <- blank_mass %>% mutate_at(c("log_masspm","day"), as.numeric)

blank_mass<-blank_mass[!(blank_mass$species=="recens" & blank_mass$host=="Alnus"),]

#match levels, because there are some sex blanks
blank_mass$sex <- factor(blank_mass$sex,
                         levels = levels(as.factor(all_synchrony_data$sex)))

####2.2.1 Main trend####
mass_pred_main<-predict(syncmod_mass_main,newdat=blank_mass,marginal=~culu +
                          idh(1+day):species + 
                          idh(1+day):host + 
                          idh(1+day):species:sex + 
                          idh(1+day):host:sex +
                          idh(1+day):species:host + 
                          idh(1+day):year:host +
                          idh(1+day):sex:species:host,
                        type="response",interval="confidence")
mass_predictions_main<-cbind(blank_mass, mass_pred_main)
#drop uneccessary rows
mass_predictions_main<-mass_predictions_main[which(mass_predictions_main$species=="antiqua" & mass_predictions_main$host=="Betula"),]
mass_predictions_main$fit <- exp(mass_predictions_main$fit)
mass_predictions_main$lwr <- exp(mass_predictions_main$lwr)
mass_predictions_main$upr <- exp(mass_predictions_main$upr)
mass_predictions_main$culu<-NULL
mass_predictions_main$log_masspm<-NULL
mass_predictions_main$species<-NULL
mass_predictions_main$host<-NULL
write.csv(mass_predictions_main,"FILE PATH/mass_predictions_main.csv", row.names = FALSE)

####2.2.2 Caterpillar species averages####
mass_pred_species<-predict(syncmod_mass_main,newdat=blank_mass,marginal=~culu +
                             #idh(1+day):species + 
                             idh(1+day):host + 
                             #idh(1+day):species:sex + 
                             idh(1+day):host:sex +
                             idh(1+day):species:host + 
                             idh(1+day):year:host +
                             idh(1+day):sex:species:host,
                           type="response",interval="confidence")
mass_predictions_species<-cbind(blank_mass, mass_pred_species)
mass_predictions_species<-mass_predictions_species[which(mass_predictions_species$host=="Betula"),]
mass_predictions_species$fit <- exp(mass_predictions_species$fit)
mass_predictions_species$lwr <- exp(mass_predictions_species$lwr)
mass_predictions_species$upr <- exp(mass_predictions_species$upr)
mass_predictions_species$culu<-NULL
mass_predictions_species$log_masspm<-NULL
mass_predictions_species$host<-NULL
write.csv(mass_predictions_species,"FILE PATH/mass_predictions_species.csv", row.names = FALSE)

####2.2.3 Host species averages####
mass_pred_hosts<-predict(syncmod_mass_main,newdat=blank_mass,marginal=~culu +
                           idh(1+day):species + 
                           #idh(1+day):host + 
                           idh(1+day):species:sex + 
                           #idh(1+day):host:sex +
                           idh(1+day):species:host + 
                           idh(1+day):year:host +
                           idh(1+day):sex:species:host,
                         type="response",interval="confidence")
mass_predictions_hosts<-cbind(blank_mass, mass_pred_hosts)
mass_predictions_hosts<-mass_predictions_hosts[which(mass_predictions_hosts$species=="antiqua"),]
mass_predictions_hosts$fit <- exp(mass_predictions_hosts$fit)
mass_predictions_hosts$lwr <- exp(mass_predictions_hosts$lwr)
mass_predictions_hosts$upr <- exp(mass_predictions_hosts$upr)
mass_predictions_hosts$culu<-NULL
mass_predictions_hosts$log_masspm<-NULL
mass_predictions_hosts$species<-NULL
write.csv(mass_predictions_hosts,"FILE PATH/mass_predictions_hosts.csv", row.names = FALSE)

####2.2.4 All species combinations####
mass_pred_all<-predict(syncmod_mass_main,newdat=blank_mass,marginal=~culu,type="response",interval="confidence")
mass_predictions_all<-cbind(blank_mass, mass_pred_all)
mass_predictions_all$fit <- exp(mass_predictions_all$fit)
mass_predictions_all$lwr <- exp(mass_predictions_all$lwr)
mass_predictions_all$upr <- exp(mass_predictions_all$upr)
mass_predictions_all$culu<-NULL
mass_predictions_all$log_masspm<-NULL
write.csv(mass_predictions_all,"FILE PATH/mass_predictions_all.csv", row.names = FALSE)


####2.3 DEVELOPMENT TIME####
#generate blank dataframe for predictions
blank_dev<-blank_mass
names(blank_dev)[names(blank_dev) == "log_masspm"] <- "durationl"

####2.3.1 Main trend####
dev_pred_main<-predict(syncmod_dev_main,newdat=blank_dev,marginal=~culu +
                         idh(1+day):species + 
                         idh(1+day):host + 
                         idh(1+day):species:sex + 
                         idh(1+day):host:sex +
                         idh(1+day):species:host + 
                         idh(1+day):year:host +
                         idh(1+day):sex:species:host,type="response",interval="confidence")
dev_predictions_main<-cbind(blank_dev, dev_pred_main)
dev_predictions_main<-dev_predictions_main[which(dev_predictions_main$species=="antiqua" & dev_predictions_main$host=="Betula"),]
dev_predictions_main$culu<-NULL
dev_predictions_main$durationl<-NULL
dev_predictions_main$species<-NULL
dev_predictions_main$host<-NULL
write.csv(dev_predictions_main,"FILE PATH/dev_predictions_main.csv", row.names = FALSE)

####2.3.2 Caterpillar species averages####
dev_pred_species<-predict(syncmod_dev_main,newdat=blank_dev,marginal=~culu +
                            #idh(1+day):species + 
                            idh(1+day):host + 
                            #idh(1+day):species:sex + 
                            idh(1+day):host:sex +
                            idh(1+day):species:host + 
                            idh(1+day):year:host +
                            idh(1+day):sex:species:host,type="response",interval="confidence")
dev_predictions_species<-cbind(blank_dev, dev_pred_species)
dev_predictions_species<-dev_predictions_species[which(dev_predictions_species$host=="Betula"),]
dev_predictions_species$culu<-NULL
dev_predictions_species$durationl<-NULL
dev_predictions_species$host<-NULL
write.csv(dev_predictions_species,"FILE PATH/dev_predictions_species.csv", row.names = FALSE)

####2.3.3 Host species averages####
dev_pred_hosts<-predict(syncmod_dev_main,newdat=blank_dev,marginal=~culu +
                          idh(1+day):species + 
                          #idh(1+day):host + 
                          idh(1+day):species:sex + 
                          #idh(1+day):host:sex +
                          idh(1+day):species:host + 
                          idh(1+day):year:host +
                          idh(1+day):sex:species:host,type="response",interval="confidence")
dev_predictions_hosts<-cbind(blank_dev, dev_pred_hosts)
dev_predictions_hosts<-dev_predictions_hosts[which(dev_predictions_hosts$species=="antiqua"),]
dev_predictions_hosts$culu<-NULL
dev_predictions_hosts$durationl<-NULL
dev_predictions_hosts$species<-NULL
write.csv(dev_predictions_hosts,"FILE PATH/dev_predictions_hosts.csv", row.names = FALSE)

####2.3.4 All species combinations####
dev_pred_all<-predict(syncmod_dev_main,newdat=blank_dev,marginal=~culu,type="response",interval="confidence")
dev_predictions_all<-cbind(blank_dev, dev_pred_all)
dev_predictions_all$culu<-NULL
dev_predictions_all$durationl<-NULL
write.csv(dev_predictions_all,"FILE PATH/dev_predictions_all.csv", row.names = FALSE)
