


hop_plot <- sam_plot_fun(model_model = sev_sam, model_raw = sev_sam_raw) +
  scale_y_discrete(labels = c("Temperature","Precipitation", "Plant biomass")) 

sevwts_ppt <- sev_sam_ppt %>%
  left_join(sev_sam_ppt_raw, by = c("lag", "variable",
                                    'dataset'))

(sevwts_pptp <- ggplot(sevwts_ppt) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    ylim(0, 1.01) +
    labs(x = "Seasons into the past", 
         y = "Cumulative seasonal weights \n posterior median and 95% BCI") +
    labs(title = "Precipitation") +
    scale_x_continuous(limits = c(1,6), 
                       breaks = seq(1, 6, by = 1),
                       labels = c(0:5)) +
    theme(axis.title.y = element_blank()))

sevwts_temp <- sev_sam_temp %>%
  left_join(sev_sam_temp_raw, by = c("lag", "variable",
                                     'dataset'))

(sevwts_tempp <- ggplot(sevwts_temp) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into the past") +
    labs(title = "Temperature") +
    scale_x_continuous(limits = c(1,6), 
                       breaks = seq(1, 6, by = 1),
                       labels = c(0:5))) +
  theme(axis.title = element_blank())

sevwts_npp <- sev_sam_npp %>%
  left_join(sev_sam_npp_raw, by = c("lag", "variable",
                                    'dataset'))

(sevwts_nppp <- ggplot(sevwts_npp) +
    geom_ribbon(aes(x = lag, ymin = LCI_obs,
                    ymax = UCI_obs),
                fill = observed_col, alpha = 0.5) +
    geom_line(aes(x = lag, y = median_obs), 
              color = observed_col, alpha = 0.5) +
    geom_pointrange(aes(x = lag, y = median,
                        ymin = LCI, 
                        ymax = UCI),
                    color = modeled_col) +    
    geom_line(aes(x = lag, y = median), 
              color = modeled_col, alpha = 0.5,
              linewidth = 0.6) +
    labs(x = "Seasons into the past", 
         y = "Cumulative seasonal \nweights") +
    labs(title = "Plant biomass") +
    scale_x_continuous(limits = c(1,11), 
                       breaks = seq(1, 11, by = 1),
                       labels = c(0:10)))

sevwts_tempp <- sevwts_tempp +
  theme(axis.title = element_blank())

sevwts_nppp <- sevwts_nppp +
  theme(axis.title.x = element_blank())
hop_plot /
(sevwts_nppp + sevwts_pptp + sevwts_tempp) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave(here('pictures',
            'sam_models',
            'grasshopper_results.pdf'
            ),
       width = 7.5,
       height = 5)
