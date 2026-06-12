####3.0 LOAD RELEVANT DATA####
all_synchrony_data <- read.csv("FILE PATH/all_synchrony_data.csv")
all_synchrony_data <- all_synchrony_data %>%
  mutate_at(vars(species,host,year,culu,sex),as.factor) %>%
  mutate_at(vars(survp,masspm,durationl,day),as.numeric)

surv_predictions_main <- read.csv("FILE PATH/surv_predictions_main.csv")
surv_predictions_hosts <- read.csv("FILE PATH/surv_predictions_hosts.csv")
surv_predictions_species <- read.csv("FILE PATH/surv_predictions_species.csv")
surv_predictions_all <- read.csv("FILE PATH/surv_predictions_all.csv")

mass_predictions_main <- read.csv("FILE PATH/mass_predictions_main.csv")
dev_predictions_main <- read.csv("FILE PATH/dev_predictions_main.csv")

####3.1 FIGURE 2####
#survival
surv_data<-all_synchrony_data[!is.na(all_synchrony_data$survp),]
surv_groups<-group_by(surv_data, species, host, day, year)
surv_summary<-summarise(surv_groups, mean(survp))
names(surv_summary)[names(surv_summary) == "mean(survp)"] <- "survp"
surv_summary<-surv_summary[which(surv_summary$year=="2021"),]

surv_subset<-as.data.frame(surv_predictions_main[which(surv_predictions_main$year=="2021"),])

surv_main_plot<-ggplot(data=surv_subset, aes(x=day, y=fit))+
  geom_point(data=surv_summary, aes(x=day, y=survp), size=1.5, stat="identity", size = 1.2,
             position=position_jitter(height = 0.04, width = 0.75), colour="palegreen3", alpha=0.25)+
  geom_ribbon(aes(x=day, ymin=lwr, ymax=upr), linetype=0, alpha=0.4, fill="palegreen4")+
  geom_line(alpha=1, size=1, colour="palegreen4")+
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(name="Asynchrony (calendar days)", limit = c(0, 65))+
  scale_y_continuous(name="Survival Probability", limit = c(0, 1, 0.25))

#mass
mass_data<-all_synchrony_data[!is.na(all_synchrony_data$masspm),]
mass_groups<-group_by(mass_data, species, host, day, sex, year)
mass_summary<-summarise(mass_groups, mean(masspm))
mass_summary<-as.data.frame(mass_summary[which(mass_summary$sex=="F" & mass_summary$year=="2021"),])
names(mass_summary)[names(mass_summary) == "mean(masspm)"] <- "masspm"

mass_subset<-as.data.frame(mass_predictions_main[which(mass_predictions_main$sex=="F" & mass_predictions_main$year=="2021"),])

mass_main_plot<-ggplot(data=mass_subset, aes(x=day, y=fit))+
  geom_point(data=mass_summary, aes(x=day, y=masspm),stat="identity", size = 1.2,
             size=1.5, position=position_jitter(height = 0, width = 0.5), colour="darkorange3", alpha=0.25)+
  geom_ribbon(aes(x=day, ymin=lwr, ymax=upr), linetype=0, alpha=0.3, fill="darkorange3")+
  geom_line(alpha=1, size=1, colour="darkorange3")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8))+
  scale_x_continuous(name="Asynchrony (days)", limit = c(0, 65))+
  scale_y_continuous(name="Pupal Mass (mg)", limit = c(0, 750, 20))

#dev
dev_data<-all_synchrony_data[!is.na(all_synchrony_data$durationl),]
dev_groups<-group_by(dev_data, species, host, day, sex, year)
dev_summary<-summarise(dev_groups, mean(durationl))
dev_summary<-as.data.frame(dev_summary[which(dev_summary$sex=="F" & dev_summary$year=="2021"),])
names(dev_summary)[names(dev_summary) == "mean(durationl)"] <- "durationl"

dev_subset<-as.data.frame(dev_predictions_main[which(dev_predictions_main$sex=="F" & dev_predictions_main$year=="2021"),])

dev_main_plot<-ggplot(data=dev_subset, aes(x=day, y=fit))+
  geom_point(data=dev_summary, aes(x=day, y=durationl),stat="identity", size = 1.2,
             position=position_jitter(height = 0, width = 0.2), colour="skyblue4", alpha=0.25)+
  geom_ribbon(aes(x=day, ymin=lwr, ymax=upr), linetype=0, alpha=0.25, fill="skyblue4")+
  geom_line(alpha=1, size=1, colour="skyblue4")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8))+
  scale_x_continuous(name="Asynchrony (days)", limit = c(0, 65))+
  scale_y_continuous(name="Dev. Time (days)", limit = c(20, 60, 10))

figure_panel_main<-ggarrange(surv_main_plot, 
                             mass_main_plot + rremove("ylab") + rremove("xlab"),
                             dev_main_plot + rremove("ylab") + rremove("xlab"),
                             labels = c("a", "b", "c"),
                             nrow=1, ncol=3,
                             align = "hv")

ggsave("FILE PATH/main_plot.svg", units = "mm",
       plot = figure_panel_main, device = "svg", width = 150, height = 45)

####3.2 FIGURE 3ab (average species plots)####
#species
surv_groups<-group_by(all_synchrony_data, day, year, species)
surv_summary<-summarise(surv_groups, mean(survp))
names(surv_summary)[names(surv_summary) == "mean(survp)"] <- "survp"
surv_summary<-surv_summary[which(surv_summary$year=="2021"),]

species_subset<-as.data.frame(surv_predictions_species[which(surv_predictions_species$year=="2021"),])

species_points1<-species_subset[which(species_subset$day=="0"),]
species_points1<- species_points1 %>%
  mutate(day = day-(2*row_number()))
species_points2<-species_subset[which(species_subset$day=="65"),]
species_points2<- species_points2 %>%
  mutate(day = day+(2*row_number()))
species_points<-rbind(species_points1,species_points2)

species_plot<-ggplot(data=species_subset, aes(x=day, y=fit, group=species, color=species))+
  annotate("rect", xmin = -13, xmax = -1, ymin = 0, ymax = 1,
           fill = "gray95", alpha = 0.5) +
  annotate("rect", xmin = 66, xmax = 77, ymin = 0, ymax = 1,
           fill = "gray95", alpha = 0.5) +
  geom_line(size=0.9)+
  geom_point(data=surv_summary, aes(x=day, y=survp, group=species, color=species), size=1.5, stat="identity", position="dodge")+
  geom_point(data=species_points, pch=15, aes(x=day, y=fit, group=species, color=species),stat="identity", position="dodge", size=1)+
  geom_errorbar(data=species_points, aes(x=day, ymin=lwr, ymax=upr), width=1, size=0.50)+
  scale_color_brewer(palette = "Dark2", name = "Caterpillars", breaks=c('brumata', 'defoliaria', 'dispar', 'monacha', 'antiqua', 'recens'), labels = c("Winter Moth", "Mottled Umber", "Gypsy Moth","Black Arches","Vapourer", "Scarce Vapourer"))+
  scale_fill_brewer(palette = "Dark2", name = "Caterpillars", breaks=c('brumata', 'defoliaria', 'dispar', 'monacha', 'antiqua', 'recens'), labels = c("Winter Moth", "Mottled Umber", "Gypsy Moth","Black Arches","Vapourer", "Scarce Vapourer"))+            
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.background=element_blank())+
  scale_x_continuous(name="Asynchrony (calendar days)", limit = c(-13, 77), breaks=seq(0, 60, 20))+
  scale_y_continuous(name="Survival Probability", limit = c(0, 1, 0.25))

#hosts
all_synchrony_data_hosts <- all_synchrony_data %>% drop_na(survp)
surv_groups<-group_by(all_synchrony_data_hosts, day, year, host)
surv_summary<-summarise(surv_groups, mean(survp))
names(surv_summary)[names(surv_summary) == "mean(survp)"] <- "survp"
surv_summary<-surv_summary[which(surv_summary$year=="2021"),]

hosts_subset<-as.data.frame(surv_predictions_hosts[which(surv_predictions_hosts$year=="2021"),])

host_points1<-hosts_subset[which(hosts_subset$day=="0"),]
host_points1<-host_points1 %>%
  mutate(day = day-(2*row_number()))
host_points2<-hosts_subset[which(hosts_subset$day=="65"),]
host_points2<-host_points2 %>%
  mutate(day = day+(2*row_number()))
host_points<-rbind(host_points1,host_points2)

hosts_plot<-ggplot(data=hosts_subset, aes(x=day, y=fit, group=host, color=host))+
  annotate("rect", xmin = -17, xmax = -1, ymin = 0, ymax = 1,
           fill = "gray95", alpha = 0.5) +
  annotate("rect", xmin = 66, xmax = 82, ymin = 0, ymax = 1,
           fill = "gray95", alpha = 0.5) +
  geom_line(size=0.9)+
  geom_point(data=surv_summary, aes(x=day, y=survp, group=host, color=host), size=1.5, stat="identity", position="dodge")+
  geom_point(data=host_points, pch=15, aes(x=day, y=fit, group=host, color=host), stat="identity", position="dodge", size=1)+
  geom_errorbar(data=host_points, aes(x=day, ymin=lwr, ymax=upr), width=1, size=0.50)+
  scale_color_d3(name = "Hosts", breaks=c('Acer', 'Alnus', 'Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Sycamore", "Alder", "Birch","Hawthorn","Apple", "Oak","Willow", "Sallow"))+
  scale_fill_d3(name = "Hosts", breaks=c('Acer', 'Alnus', 'Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Sycamore", "Alder", "Birch","Hawthorn","Apple", "Oak","Willow", "Sallow"))+            
  theme_bw()+
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(name="Asynchrony (calendar days)", limit = c(-17, 82), breaks=seq(0, 60, 20))+
  scale_y_continuous(name="Survival Probability", limit = c(0, 1, 0.25))

hosts_and_species_overall<-ggarrange(species_plot + rremove("ylab"), 
                                     hosts_plot,
                                     labels = c("a", "b"),
                                     nrow=2, ncol=1,
                                     legend = "right",
                                     align = "hv")

ggsave("FILE PATH/hosts_and_species_overall.svg", units = "mm",
       plot = hosts_and_species_overall, device = "svg", width = 90, height = 100)

####3.3 FIGURE 3c-h (species-specific plots)####
#validate with datapoints of raw data on graph
surv_groups<-group_by(all_synchrony_data, host, species, day)
surv_summary<-summarise(surv_groups, mean(survp))
#remove Prunus#
surv_summary<-surv_summary[!(surv_summary$host=="Prunus"),]
names(surv_summary)[names(surv_summary) == "mean(survp)"] <- "survp"

surv_summary_specific<-as.data.frame(surv_summary[which(surv_summary$species=="recens"),])

#plot
surv_predictions_subset<-surv_predictions_all[which(surv_predictions_all$species=="recens" & surv_predictions_all$year=="2021"),]

surv_recens_plot<-ggplot(data=surv_predictions_subset, aes(x=day, y=fit, group=host, color=host, fill=host))+
  geom_ribbon(aes(x=day, ymin=lwr, ymax=upr), linetype=0, alpha=0.25)+
  geom_line(alpha=1, size=0.9)+
  geom_point(data=surv_summary_specific, aes(x=day, y=survp, group=host, color=host), size=1.5, stat="identity", position="dodge")+
  scale_x_continuous(name="Asynchrony (days)", limit = c(0, 65))+
  scale_y_continuous(name="Survival Probability", limit = c(0, 1, 0.25))+
  theme_bw()+
  scale_color_d3(name = "Host-plant Species", breaks=c('Acer', 'Alnus','Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Acer", "Alnus", "Betula", "Crataegus", "Malus", "Quercus", "S. alba", "S. caprea"))+
  scale_fill_d3(name = "Host-plant Species", breaks=c('Acer', 'Alnus','Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Acer", "Alnus", "Betula", "Crataegus", "Malus", "Quercus", "S. alba", "S. caprea"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8))

#annotations as x=22 and y=0.1 for O. antiqua#
#antiqua x=28, y=0.1
#dispar x=60, y=0.95
#recens x=60, y=0.95
#monacha x=60, y=0.95
#brumata x=60, y=0.95
#defoliaria x=60, y=0.95

#plotting#
surv_host_by_species<-ggarrange(surv_brumata_plot + rremove("ylab") + rremove("xlab"), 
                                surv_defoliaria_plot + rremove("ylab") + rremove("xlab"),
                                surv_dispar_plot + rremove("ylab") + rremove("xlab"),
                                surv_monacha_plot + rremove("ylab") + rremove("xlab"),
                                surv_antiqua_plot,
                                surv_recens_plot + rremove("ylab") + rremove("xlab"),
                                labels = c("c", "d", "e", "f", "g", "h"),
                                nrow=3, ncol=2,
                                common.legend = TRUE, legend = "none",
                                align = "hv")

ggsave("FILE PLOT/surv_host_by_species.svg", units = "mm",
       plot = surv_host_by_species, device = "svg", width = 90, height = 130)