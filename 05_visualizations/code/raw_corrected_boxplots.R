#observed vs. corrected dissimilarity analysis
#Ana Miller-ter Kuile
#October 16, 2023

#this is an exploratory script to get the paired differences
#among observed and corrected dissimilarlity metrics, but 
#looking at differences among datasets too

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  'emmeans', 'glmmTMB',
                  'patchwork')


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Simulated process -------------------------------------------------------
# 
# type <- c(rep("observed", 10), rep("corrected", 10))
# samp <- c(1:10, 1:10)
# dataset <- c(rep(1,3), rep(2, 3), rep(3,4),rep(1,3), rep(2, 3), rep(3,4))
# response <- c(runif(10, min = 0, max = 1), runif(10, min = 0, max = 0.7))
# 
# df <- as.data.frame(cbind(type = type,
#                     samp = samp,
#                     dataset = dataset,
#                     response = response)) %>%
#   mutate(response = as.numeric(response))
# 
# 
# m <- glmmTMB(response ~ type*dataset + (1|samp),
#              data = df,
#              beta_family)
# 
# summary(m)
# 
# ggplot(df, aes(x = dataset, y = response, fill = type)) +
#   geom_boxplot()

#could also run this as a JAGS model - but for now this works I think
#as a template for assessing this.


# Load data ---------------------------------------------------------------

sbc_obs <- readRDS(here('05_visualizations',
                        'viz_data',
                        'sbc_observed_bray.RDS'))

str(sbc_obs)

sbc_modeled <- readRDS(here('01_sbc_fish',
                            'monsoon',
                            'fish_MSAM',
                            'outputs',
                            'fish_bray_meanSD.RDS'))
#will need:
#raw bray for all communities, 
#corrected bray for communities - linked to raw

konza_obs <- readRDS(here('05_visualizations',
                          'viz_data',
                          'konza_observed_bray.RDS'))

konza_modeled <- readRDS(here('02_konza_birds',
                              'monsoon',
                              'MSAM',
                              'outputs',
                              'bird_bray_meanSD.RDS'))

#sevilleta
sev_obs <- readRDS(here('05_visualizations',
                        'viz_data',
                        'sev_observed_bray.RDS'))

sev_modeled <- readRDS(here('03_sev_grasshoppers',
                            'monsoon',
                            "MSAM",
                            "outputs",
                            'sev_bray_meanSD.RDS'))

nps_obs <- readRDS(here('05_visualizations',
                        'viz_data',
                        'nps_observed_jaccard.RDS'))

nps_obs <- nps_obs %>%
  rename(yrID = EventYear,
         siteID = plot_trans_quad)

nps_modeled <- readRDS(here('04_nps_plants',
                            'monsoon',
                            'nps_MSAM',
                            'outputs_yrsite',
                            'nps_Jaccard_summary.RDS'))

# Prep modeled data -------------------------------------------------------

#site x year
sbc_m2 <- as.data.frame(sbc_modeled) %>%
  rownames_to_column(var = "var") %>%
  filter(var != "deviance") %>%
  separate(var,
           into = c('siteID', 'yrID'),
           sep = ",") %>%
  mutate(siteID = str_sub(siteID, 6, nchar(siteID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  rename("bray" = "Mean") %>%
  dplyr::select(yrID, siteID, bray) %>%
  mutate(type = "modeled")%>%
  mutate(yrID = as.numeric(yrID),
         siteID = as.numeric(siteID))

str(sbc_m2)

kz_m2 <- as.data.frame(konza_modeled) %>%
  rownames_to_column(var = "var") %>%
  filter(var != "deviance") %>%
  separate(var,
           into = c('siteID', 'yrID'),
           sep = ",") %>%
  mutate(siteID = str_sub(siteID, 6, nchar(siteID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  rename("bray" = "Mean") %>%
  dplyr::select(yrID, siteID, bray) %>%
  mutate(type = "modeled")%>%
  mutate(yrID = as.numeric(yrID),
         siteID = as.numeric(siteID))
  
str(kz_m2)

sev_m2 <- as.data.frame(sev_modeled) %>%
  rownames_to_column(var = "var") %>%
  filter(var != "deviance") %>%
  separate(var, 
           into = c("siteID", "yrID"),
           sep = ",")%>%
  mutate(siteID = str_sub(siteID, 6, nchar(siteID))) %>%
  mutate(yrID = str_sub(yrID, 1, (nchar(yrID)-1))) %>%
  rename("bray" = "Mean") %>%
  dplyr::select(yrID, siteID, bray) %>%
  mutate(type = "modeled")%>%
  mutate(yrID = as.numeric(yrID),
         siteID = as.numeric(siteID))

str(sev_m2)
nps_m2 <- nps_modeled %>%
  rename(yrID = EventYear,
         siteID = plot_trans_quad,
         turnover = mean) %>%
  dplyr::select(-sd, - mean_loss, -sd_loss, -sd_gain, -mean_gain) %>%
  mutate(type = "modeled")

str(nps_m2)

# Combine -----------------------------------------------------------------

sbc_bray <- rbind(sbc_obs, sbc_m2) %>%
  filter(!is.na(bray)) %>%
  unite("site_year",
        c(siteID, yrID),
        sep = "_",
        remove = F) %>%
  #the beta family in glmmTMB doesn't work
  #if values are exactly 1 or exactly 0
  mutate(bray = case_when(bray == 0 ~ 0.001,
                          bray == 1 ~ 0.9999,
                          TRUE ~ bray)) %>%
  mutate(dataset = "sbc_fish") %>%
  rename(diss = bray)

kz_bray <- konza_obs %>%
  rename("siteID" = "TransID") %>% 
  rbind(kz_m2) %>%
  filter(!is.na(bray)) %>%
  unite("site_year",
        c(siteID, yrID),
        sep = "_",
        remove = F) %>%
  #the beta family in glmmTMB doesn't work
  #if values are exactly 1 or exactly 0
  mutate(bray = case_when(bray == 0 ~ 0.001,
                          bray == 1 ~ 0.9999,
                          TRUE ~ bray)) %>%
  mutate(dataset = "konza_birds") %>%
  rename(diss = bray)

sev_bray <- sev_obs %>%
  rbind(sev_m2)  %>%
  unite("site_year",
        c(siteID, yrID),
        sep = "_",
        remove = F) %>%
  #the beta family in glmmTMB doesn't work
  #if values are exactly 1 or exactly 0
  mutate(bray = case_when(bray == 0 ~ 0.001,
                          bray == 1 ~ 0.9999,
                          TRUE ~ bray)) %>%
  mutate(dataset = "sev_hoppers") %>%
  rename(diss = bray)

nps_turn <- nps_obs %>%
  rbind(nps_m2)%>% 
  filter(!is.na(turnover)) %>%
  unite("site_year",
        c(siteID, yrID),
        sep = "_",
        remove = F) %>%
  #the beta family in glmmTMB doesn't work
  #if values are exactly 1 or exactly 0
  mutate(turnover = case_when(turnover == 0 ~ 0.001,
                          turnover == 1 ~ 0.9999,
                          TRUE ~ turnover)) %>%
  mutate(dataset = "nps_plants") %>%
  rename(diss = turnover)


all_diss <- rbind(sbc_bray, kz_bray, sev_bray, nps_turn) %>% 
  # making type labels nice
  mutate(type_label = case_match(type, 
    "modeled" ~ "Modeled\n(imperfect detection)",
    "observed_all" ~ "Empirical\n(maximum)",
    "observed_one" ~ "Empirical\n(single survey)"
  ))


# Visualize ---------------------------------------------------------------


modeled_col <- "#E88C23"
observed_col <- "#438AA8"

boxplot_function <- function(dataset) {
  
  if(dataset == "birds"){
    df <- all_diss %>% 
      filter(dataset == "konza_birds")
  } else if(dataset == "fish") {
    df <- all_diss %>% 
      filter(dataset == "sbc_fish")
  } else if(dataset == "grasshoppers"){
    df <- all_diss %>% 
      filter(dataset == "sev_hoppers")
  } else if(dataset == "plants"){
    df <- all_diss %>% 
      filter(dataset == "nps_plants")
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    title = "(b) Birds"
  } else if(dataset == "fish") {
    title = "(a) Fish"
  } else if(dataset == "grasshoppers") {
    title = "(c) Grasshoppers"
  } else if(dataset == "plants") {
    title = "(d) Plants"
  }else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  df %>% 
    ggplot(aes(x = type_label, y = 1-diss, fill = type)) +
    geom_violin(aes(color = type, fill = type), alpha = 0.7) +
    scale_fill_manual(values = c(modeled = modeled_col, 
                                 observed_all = observed_col, 
                                 observed_one = observed_col)) +
    scale_color_manual(values = c(modeled = modeled_col, 
                                 observed_all = observed_col, 
                                 observed_one = observed_col)) +
    geom_boxplot(aes(fill = type), width = 0.1, outlier.size = 1, alpha = 0.9) +
    labs(x = "Type", y = "Stability\n(1 - dissimilarity metric)", title = title) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none",
          plot.title.position = "plot",
          panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8)) 
  
}

knz_boxplot <- boxplot_function("birds") 
knz_boxplot

sbc_boxplot <- boxplot_function("fish") 
sbc_boxplot

sev_boxplot <- boxplot_function("grasshoppers")
sev_boxplot

nps_boxplot <- boxplot_function("plants")
nps_boxplot

all_boxplot <- (sbc_boxplot | knz_boxplot) /
               (sev_boxplot | nps_boxplot)
all_boxplot 


ggsave(plot = all_boxplot,
       filename = here("pictures",
                       "detection_models",
                       "observed_modeled_violin.jpg"),
       height = 12,
       width = 16,
       units = "cm",
       dpi = 400)

# Looking at differences across datasets ----------------------------------

library(brms)


# a1 <- brm(log_yield ~ log_width + NHD_RdDensWs +
#             Dam_binary + meanTemp + exc_y + (1|huc_2), 
#           data = dat_amax_brms, family = gaussian())


m1 <- brm(diss ~ type*dataset + (1|site_year),
          data = all_diss,
          family = Beta())

summary(m1)

em <- emmeans(m1, pairwise ~ type | dataset)

em

t <- as.data.frame(em$contrasts)

ggplot(t, aes(x = dataset, y= estimate, shape = contrast)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), 
                width = 0, linewidth = 0.75,
                position = position_dodge(width = 0.4))

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "detection_models",
                       "observed_modeled_data_contrasts.jpg"),
       height = 4,
       width = 6,
       units = "in")
