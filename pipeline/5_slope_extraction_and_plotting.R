####5.0 EXTRACTING SLOPES AND STACKING THEM####

get_day_effects <- function(model, species_vec, host_vec, year_vec) {
  
  # Extract Sol matrix
  sol <- model$Sol
  
  # Empty list to store results
  results <- list()
  
  # Loop through all combinations
  for (sp in species_vec) {
    for (ho in host_vec) {
      for (yr in year_vec) {
        
        # Build the term names according to Sol column names
        terms <- c("day")  # main effect always included
        
        sp_term <- paste0("day.species.", sp)
        if (sp_term %in% colnames(sol)) terms <- c(terms, sp_term)
        
        ho_term <- paste0("day.host.", ho)
        if (ho_term %in% colnames(sol)) terms <- c(terms, ho_term)
        
        sp_ho_term <- paste0("day.species:host.", sp, ".", ho)
        if (sp_ho_term %in% colnames(sol)) terms <- c(terms, sp_ho_term)
        
        yr_term <- paste0("day:year", yr)
        if (yr_term %in% colnames(sol)) terms <- c(terms, yr_term)
        
        yr_ho_term <- paste0("day.year:host.", yr, ".", ho)
        if (yr_ho_term %in% colnames(sol)) terms <- c(terms, yr_ho_term)
        
        # Sum posterior samples for this combination
        combined_effect <- rowSums(sol[, terms, drop=FALSE])
        
        # Calculate mean and HPD interval
        mean_val <- mean(combined_effect)
        hpd <- HPDinterval(as.mcmc(combined_effect))
        
        # Store results
        results[[paste(sp, ho, yr, sep="_")]] <- data.frame(
          species = sp,
          host = ho,
          year = yr,
          mean = mean_val,
          lower = hpd[1],
          upper = hpd[2]
        )
      }
    }
  }
  
  # Combine into a single dataframe
  do.call(rbind, results)
}

species_vec <- c("brumata", "defoliaria", "dispar", "monacha", "antiqua", "recens")
host_vec <- c("Acer", "Alnus", "Betula", "Crataegus", "Malus", "Quercus", "Salba", "Scapraea")
year_vec <- c("2021")

results_df <- get_day_effects(syncmod_surv_main, species_vec, host_vec, year_vec)

rownames(results_df) <- NULL

write.csv(results_df, file="FILE PATH/surv_slopes.csv", row.names=FALSE)
results_df <- read.csv("FILE PATH/surv_slopes.csv")

#arrange low to high

results_df <- results_df %>%
  arrange(mean) %>%
  mutate(order = row_number())

#plot slopes to show which species are steepest and which are shallowest

results_df$species <- factor(
  results_df$species,
  levels = c("brumata", "defoliaria", "dispar", "monacha", "antiqua", "recens")
)

surv_slope_plot<-ggplot(data=results_df, aes(x=species, y=mean, color=host))+
  geom_point(position=position_dodge(0.5))+
  geom_errorbar(aes(ymin=lower, ymax=upper, colour=host), position=position_dodge(0.5), width=0)+
  theme_bw()+
  scale_color_d3(name = "Hosts", breaks=c('Acer', 'Alnus', 'Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Sycamore", "Alder", "Birch","Hawthorn","Apple", "Oak","Willow", "Sallow"))+
  scale_fill_d3(name = "Hosts", breaks=c('Acer', 'Alnus', 'Betula', 'Crataegus', 'Malus', 'Quercus', 'Salba', 'Scapraea'), labels = c("Sycamore", "Alder", "Birch","Hawthorn","Apple", "Oak","Willow", "Sallow"))+            
  theme(axis.title.x = element_text(face="bold"), axis.title.y = element_text(face="bold"), 
        legend.title=element_text(face="bold"), text = element_text(family = "Sans", size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))+
  scale_x_discrete(name="Caterpillar Species", labels = c("Winter Moth", "Mottled Umber", "Gypsy Moth", "Black Arches", "Vapourer", "Scarce Vapourer"))+
  scale_y_continuous(name="Survival Slope (log odds/day)", limit = c(-0.13,0.01))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")

ggsave("FILE PATH/surv_slope_plot.svg", units = "mm",
       plot = surv_slope_plot, device = "svg", width = 120, height = 60)