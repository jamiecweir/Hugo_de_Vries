#The code below estimates the variance in the response (on link scale) associated with
#fixed effects, random slopes, and random intercepts in the MAIN MODEL (see 1_modelling.R)

library(MCMCglmm)

sync_day<-read.csv("PATH/all_synchrony_data.txt",sep="\t")

#### 4.1 SURVIVAL MODEL####

model<-syncmod_surv_main

#check number of fixed effects first
fixed_effects<-4

####4.1.1 Variance attributable to fixed effects####
Vx<-cov(as.matrix(model $X[,2:fixed_effects]))

Vfixed_effects<-apply(model $Sol[,2: fixed_effects], 1, function(x){t(x)%*%Vx%*%x})
median(Vfixed_effects)

#variance attribtable to day alone in 2020
median(model $Sol[,"day"]^2)*Vx["day","day"]

#variance attribtable to day deviation in 2021
median((model $Sol[,"day:year2021"]^2)*Vx["day:year2021","day:year2021"])

####4.1.2 Variance attributable to random slopes####
mean_day<-mean(sync_day$day)^2
dayslopevariance_species<-model $VCV[,"day.species"]*(Vx["day","day"]+ mean_day)
dayslopevariance_host<-model $VCV[,"day.host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_specieshost<-model $VCV[,"day.species:host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_yearhost<-model $VCV[,"day.year:host"]*(Vx["day","day"]+ mean_day)
totalslopevar<-dayslopevariance_species+ dayslopevariance_host+ dayslopevariance_specieshost+ dayslopevariance_yearhost

####4.1.3 Variance attributable to random intercepts####
totalintervar<-model$VCV[,"culu"]+model$VCV[,"(Intercept).species"]+model$VCV[,"(Intercept).host"]+model$VCV[,"(Intercept).species:host"]+model$VCV[,"(Intercept).year:host"]+model$VCV[,"units"]

####4.1.4 Total variance####
totalmodelvar<-Vfixed_effects+ totalslopevar+ totalintervar

####4.1.5 Generate variance comp dataset for survival####
variance.head<-c("Total variance","Total variance lo","Total variance hi",
                 "Fixed effect prop","Fixed effect prop lo","Fixed effect prop hi",
                 "Random intercept prop","Random intercept prop lo","Random intercept prop hi",
                 "Random slope prop","Random slope prop lo","Random slope prop hi",
                 "Random intercept species prop","Random intercept species prop lo","Random intercept species prop hi",
                 "Random intercept host prop","Random intercept host prop lo","Random intercept host prop hi",
                 "Random intercept interaction prop","Random intercept interaction prop lo","Random intercept interaction prop hi",
                 "Random slope species prop","Random slope species prop lo","Random slope species prop hi",
                 "Random slope host prop","Random slope host prop lo","Random slope host prop hi",
                 "Random slope interaction prop","Random slope interaction prop lo","Random slope interaction prop hi")

survival_outputs<-c(median(totalmodelvar),HPDinterval(totalmodelvar),
                    median(Vfixed_effects/totalmodelvar),HPDinterval(Vfixed_effects/totalmodelvar),
                    median(totalintervar/totalmodelvar),HPDinterval(totalintervar/totalmodelvar),
                    median(totalslopevar/totalmodelvar),HPDinterval(totalslopevar/totalmodelvar),
                    median(model$VCV[,"(Intercept).species"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species"]/totalintervar),
                    median(model$VCV[,"(Intercept).host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).host"]/totalintervar),
                    median(model$VCV[,"(Intercept).species:host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species:host"]/totalintervar),
                    median(dayslopevariance_species/totalslopevar),HPDinterval(dayslopevariance_species/totalslopevar),
                    median(dayslopevariance_host/totalslopevar),HPDinterval(dayslopevariance_host/totalslopevar),
                    median(dayslopevariance_specieshost/totalslopevar),HPDinterval(dayslopevariance_specieshost/totalslopevar))

####4.2 MASS MODEL####
model<-syncmod_mass_main

#check number of fixed effects first
fixed_effects<-6

####4.2.1 Variance attributable to fixed effects####
Vx<-cov(as.matrix(model $X[,2:fixed_effects]))

Vfixed_effects<-apply(model $Sol[,2: fixed_effects], 1, function(x){t(x)%*%Vx%*%x})

#variance attribtable to day alone in 2020
median(model $Sol[,"day"]^2)*Vx["day","day"]

#variance attribtable to day deviation in 2021
median((model $Sol[,"day:year2021"]^2)*Vx["day:year2021","day:year2021"])

####4.2.2 Variance attributable to random slopes####
mean_day<-mean(sync_day$day)^2

dayslopevariance_species<-model $VCV[,"day.species"]*(Vx["day","day"]+ mean_day)
dayslopevariance_host<-model $VCV[,"day.host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_specieshost<-model $VCV[,"day.species:host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_yearhost<-model $VCV[,"day.year:host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_speciessex<-model $VCV[,"day.species:sex"]*(Vx["day","day"]+ mean_day)
dayslopevariance_hostsex<-model $VCV[,"day.host:sex"]*(Vx["day","day"]+ mean_day)
dayslopevariance_speciessexhost<-model $VCV[,"day.sex:species:host"]*(Vx["day","day"]+ mean_day)


totalslopevar<-dayslopevariance_species+ dayslopevariance_host+ dayslopevariance_specieshost+ dayslopevariance_yearhost+ dayslopevariance_speciessex+ dayslopevariance_hostsex+ dayslopevariance_speciessexhost

####4.2.3 Variance attributable to random intercepts####
totalintervar<-model$VCV[,"culu"]+model$VCV[,"(Intercept).species"]+model$VCV[,"(Intercept).host"]+model$VCV[,"(Intercept).species:host"]+model$VCV[,"(Intercept).year:host"]+model$VCV[,"(Intercept).species:sex"]+model$VCV[,"(Intercept).host:sex"]+model$VCV[,"(Intercept).sex:species:host"]+model$VCV[,"units"]

####4.2.4 Total variance####
totalmodelvar<-Vfixed_effects+ totalslopevar+ totalintervar

####4.2.5 Generate variance comp dataset for mass####
mass_outputs<-c(median(totalmodelvar),HPDinterval(totalmodelvar),
                median(Vfixed_effects/totalmodelvar),HPDinterval(Vfixed_effects/totalmodelvar),
                median(totalintervar/totalmodelvar),HPDinterval(totalintervar/totalmodelvar),
                median(totalslopevar/totalmodelvar),HPDinterval(totalslopevar/totalmodelvar),
                median(model$VCV[,"(Intercept).species"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species"]/totalintervar),
                median(model$VCV[,"(Intercept).host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).host"]/totalintervar),
                median(model$VCV[,"(Intercept).species:host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species:host"]/totalintervar),
                median(dayslopevariance_species/totalslopevar),HPDinterval(dayslopevariance_species/totalslopevar),
                median(dayslopevariance_host/totalslopevar),HPDinterval(dayslopevariance_host/totalslopevar),
                median(dayslopevariance_specieshost/totalslopevar),HPDinterval(dayslopevariance_specieshost/totalslopevar))

####4.3 DEVELOPMENT MODEL####
model<-syncmod_dev_main

#check number of fixed effects first
fixed_effects<-6
####4.3.1 Variance attributable to fixed effects####
Vx<-cov(as.matrix(model $X[,2:fixed_effects]))

Vfixed_effects<-apply(model $Sol[,2: fixed_effects], 1, function(x){t(x)%*%Vx%*%x})

#variance attribtable to day alone in 2020
median(model $Sol[,"day"]^2)*Vx["day","day"]

#variance attribtable to day deviation in 2021
median((model $Sol[,"day:year2021"]^2)*Vx["day:year2021","day:year2021"])

####4.3.2 Variance attributable to random slopes####
mean_day<-mean(sync_day$day)^2

dayslopevariance_species<-model $VCV[,"day.species"]*(Vx["day","day"]+ mean_day)
dayslopevariance_host<-model $VCV[,"day.host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_specieshost<-model $VCV[,"day.species:host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_yearhost<-model $VCV[,"day.year:host"]*(Vx["day","day"]+ mean_day)
dayslopevariance_speciessex<-model $VCV[,"day.species:sex"]*(Vx["day","day"]+ mean_day)
dayslopevariance_hostsex<-model $VCV[,"day.host:sex"]*(Vx["day","day"]+ mean_day)
dayslopevariance_speciessexhost<-model $VCV[,"day.sex:species:host"]*(Vx["day","day"]+ mean_day)


totalslopevar<-dayslopevariance_species+ dayslopevariance_host+ dayslopevariance_specieshost+ dayslopevariance_yearhost+ dayslopevariance_speciessex+ dayslopevariance_hostsex+ dayslopevariance_speciessexhost

####4.3.3 Variance attributable to random intercepts####
totalintervar<-model$VCV[,"culu"]+model$VCV[,"(Intercept).species"]+model$VCV[,"(Intercept).host"]+model$VCV[,"(Intercept).species:host"]+model$VCV[,"(Intercept).year:host"]+model$VCV[,"(Intercept).species:sex"]+model$VCV[,"(Intercept).host:sex"]+model$VCV[,"(Intercept).sex:species:host"]+model$VCV[,"units"]

####4.3.4 Total variance####
totalmodelvar<-Vfixed_effects+ totalslopevar+ totalintervar

####4.3.5 Generate variance comp dataset for mass####
dev_outputs<-c(median(totalmodelvar),HPDinterval(totalmodelvar),
               median(Vfixed_effects/totalmodelvar),HPDinterval(Vfixed_effects/totalmodelvar),
               median(totalintervar/totalmodelvar),HPDinterval(totalintervar/totalmodelvar),
               median(totalslopevar/totalmodelvar),HPDinterval(totalslopevar/totalmodelvar),
               median(model$VCV[,"(Intercept).species"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species"]/totalintervar),
               median(model$VCV[,"(Intercept).host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).host"]/totalintervar),
               median(model$VCV[,"(Intercept).species:host"]/totalintervar),HPDinterval(model$VCV[,"(Intercept).species:host"]/totalintervar),
               median(dayslopevariance_species/totalslopevar),HPDinterval(dayslopevariance_species/totalslopevar),
               median(dayslopevariance_host/totalslopevar),HPDinterval(dayslopevariance_host/totalslopevar),
               median(dayslopevariance_specieshost/totalslopevar),HPDinterval(dayslopevariance_specieshost/totalslopevar))

output<-as.data.frame(rbind(survival_outputs,mass_outputs,dev_outputs))
names(output)<-variance.head

row.names(output)<-c("survival","mass","duration")

write.table(output,"FILE PATH/variance.decomp.txt",sep="\t",col.names=TRUE,row.names=TRUE)

#### 4.4 PLOTTING####

variance_decomp_tabulated <- read.csv("FILE PATH/variance_decomp_tabulated.csv")

variance_decomp_tabulated <- variance_decomp_tabulated %>%
  mutate_at(vars(model,term),as.factor) %>%
  mutate_at(vars(mean,lwr,upr),as.numeric)

variance_decomp_tabulated <- variance_decomp_tabulated %>%
  mutate(term_type = case_when(
    term %in% c("fixed", "random int", "random slope") ~ "overall",
    term %in% c("sp int", "hp int", "int int") ~ "int",
    term %in% c("sp slope", "hp slope", "int slope") ~ "slope",
    term %in% c("total") ~ "total",
  ))

variance_decomp_tabulated$term <- factor(variance_decomp_tabulated$term, levels = c("total", "fixed", "random int", "random slope", "sp int",
                                                                                    "hp int", "int int", "sp slope", "hp slope", "int slope"))

variance_decomp_tabulated$term_type <- factor(variance_decomp_tabulated$term_type, levels = c("overall", "int", "slope"))
variance_decomp_tabulated<-variance_decomp_tabulated[!(variance_decomp_tabulated$term %in% "total"),]

surv_var_comp<-variance_decomp_tabulated[which(variance_decomp_tabulated$model=="survival"),]
mass_var_comp<-variance_decomp_tabulated[which(variance_decomp_tabulated$model=="mass"),]
dev_var_comp<-variance_decomp_tabulated[which(variance_decomp_tabulated$model=="duration"),]

surv_var_plot<-ggplot(surv_var_comp, aes(x=term, y=mean, color=term_type))+
  geom_point(position=position_dodge(0.5), size=2)+
  geom_errorbar(aes(x=term, ymin=lwr, ymax=upr, color=term_type), position=position_dodge(0.5), width=0)+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        llegend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_brewer(palette = "Dark2", name = NULL, labels=c("overall"="Overall",
                                                              "int"="Intercept variance",
                                                              "slope"="Slope variance"))+
  scale_x_discrete(name= "Term", 
                   guide = guide_axis(angle = 45),
                   labels=c("fixed" = "Fixed Effects", "random int" = "Random Intercepts",
                            "random slope" = "Random Slopes", "sp int" = "Caterpillar",
                            "hp int" = "Host-plant", "int int" = "Interaction",
                            "sp slope" = "Caterpillar", "hp slope" = "Host-plant",
                            "int slope" = "Interaction"))+
  scale_y_continuous(name="Variance Explained (%)", guide = guide_axis(angle = 45))+
  coord_flip()+
  ggtitle("Survival")

mass_var_plot<-ggplot(mass_var_comp, aes(x=term, y=mean, color=term_type))+
  geom_point(position=position_dodge(0.5), size=2)+
  geom_errorbar(aes(x=term, ymin=lwr, ymax=upr, color=term_type), position=position_dodge(0.5), width=0)+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_brewer(palette = "Dark2", name = NULL, labels=c("overall"="Overall",
                                                              "int"="Intercept variance",
                                                              "slope"="Slope variance"))+
  scale_x_discrete(name= "Term", 
                   guide = guide_axis(angle = 45),
                   labels=c("fixed" = "Fixed Effects", "random int" = "Random Intercepts",
                            "random slope" = "Random Slopes", "sp int" = "Caterpillar",
                            "hp int" = "Host-plant", "int int" = "Interaction",
                            "sp slope" = "Caterpillar", "hp slope" = "Host-plant",
                            "int slope" = "Interaction"))+
  scale_y_continuous(name="Variance Explained (%)", guide = guide_axis(angle = 45))+
  coord_flip()+
  ggtitle("Mass")

dev_var_plot<-ggplot(dev_var_comp, aes(x=term, y=mean, color=term_type))+
  geom_point(position=position_dodge(0.5), size=2)+
  geom_errorbar(aes(x=term, ymin=lwr, ymax=upr, color=term_type), position=position_dodge(0.5), width=0)+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_brewer(palette = "Dark2", name = NULL, labels=c("overall"="Overall",
                                                              "int"="Intercept variance",
                                                              "slope"="Slope variance"))+
  scale_x_discrete(name= "Term", 
                   guide = guide_axis(angle = 45),
                   labels=c("fixed" = "Fixed Effects", "random int" = "Random Intercepts",
                            "random slope" = "Random Slopes", "sp int" = "Caterpillar",
                            "hp int" = "Host-plant", "int int" = "Interaction",
                            "sp slope" = "Caterpillar", "hp slope" = "Host-plant",
                            "int slope" = "Interaction"))+
  scale_y_continuous(name="Variance Explained (%)", guide = guide_axis(angle = 45))+
  coord_flip()+
  ggtitle("Duration")

var_comp<-ggarrange(surv_var_plot + rremove("xlab") + rremove("ylab"), 
                    mass_var_plot + rremove("xlab") + rremove("ylab"),
                    dev_var_plot,
                    labels = c("a", "b", "c"),
                    nrow=3, ncol=1,
                    common.legend = TRUE, legend = "right",
                    align = "hv")

ggsave("FILE PATH/var_comp_final.svg", 
       plot=var_comp, device="svg", width = 85, height = 135, units = "mm")

