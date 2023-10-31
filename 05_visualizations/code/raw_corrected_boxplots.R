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

type <- c(rep("observed", 10), rep("corrected", 10))
samp <- c(1:10, 1:10)
dataset <- c(rep(1,3), rep(2, 3), rep(3,4),rep(1,3), rep(2, 3), rep(3,4))
response <- c(runif(10, min = 0, max = 1), runif(10, min = 0, max = 0.7))

df <- as.data.frame(cbind(type = type,
                    samp = samp,
                    dataset = dataset,
                    response = response)) %>%
  mutate(response = as.numeric(response))


m <- glmmTMB(response ~ type*dataset + (1|samp),
             data = df,
             beta_family)

summary(m)

ggplot(df, aes(x = dataset, y = response, fill = type)) +
  geom_boxplot()

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
  mutate(dataset = "sbc_fish")

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
  mutate(dataset = "konza_birds")

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
  mutate(dataset = "sev_hoppers")

all_bray <- rbind(sbc_bray, kz_bray, sev_bray)


# Visualize ---------------------------------------------------------------


modeled_col <- "#E88C23"
observed_col <- "#438AA8"

ggplot(all_bray, aes(x = dataset, y = bray)) +
  #geom_jitter(aes(group = type, color = type), alpha = 0.2, width = 0.2) +
  geom_boxplot(aes(color = type),  outlier.shape = NA) + 
  geom_point(aes(color = type), position = position_jitterdodge())+ 
  labs(x = "Dataset",
       y = "Dissimilarity") +
  annotate(geom = "text", x = 0.75, y = 1, label = "More different") +
  annotate(geom = "text", x = 0.75, y = 0, label = "More similar") +
  scale_fill_manual(values = c(NA, NA)) + 
  scale_color_manual(values = c(modeled = modeled_col, observed = observed_col))  

ggplot(all_bray, aes(x = type, y = bray)) +
  #geom_jitter(aes(group = type, color = type), alpha = 0.2, width = 0.2) +
  geom_violin(aes(fill = type)) + 
 # geom_point(aes(color = type), position = position_jitterdodge())+ 
  labs(x = "Dataset",
       y = "Dissimilarity") +
  annotate(geom = "text", x = 0.75, y = 1, label = "More different") +
  annotate(geom = "text", x = 0.75, y = 0, label = "More similar") +
  scale_fill_manual(values = c(modeled = modeled_col, observed = observed_col)) +
  facet_wrap(~ dataset)

boxplot_function <- function(dataset) {
  
  if(dataset == "birds"){
    df <- all_bray %>% 
      filter(dataset == "konza_birds")
  } else if(dataset == "fish") {
    df <- all_bray %>% 
      filter(dataset == "sbc_fish")
  } else if(dataset == "grasshoppers"){
    df <- all_bray %>% 
      filter(dataset == "sev_hoppers")
  } else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  if(dataset == "birds"){
    title = "KNZ birds"
  } else if(dataset == "fish") {
    title = "SBC fish"
  } else if(dataset == "grasshoppers") {
    title = "SEV grasshoppers"
  }else {
    warning("Check your arguments! You may have specified the wrong dataset.")
    return(NA)
  }
  
  df %>% 
    ggplot(aes(x = type, y = bray, fill = type)) +
    geom_violin() +
    scale_fill_manual(values = c(modeled = modeled_col, observed = observed_col)) +
    geom_boxplot(width = 0.1) +
    labs(x = "Type", y = "Bray-Curtis dissimilarity", title = title) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none",
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)) 
  
}

knz_boxplot <- boxplot_function("birds") +
  annotate(geom = "text", x = 0.75, y = 1, label = "More different") +
  annotate(geom = "text", x = 0.75, y = 0, label = "More similar") 
knz_boxplot

sbc_boxplot <- boxplot_function("fish")
sbc_boxplot

sev_boxplot <- boxplot_function("grasshoppers")
sev_boxplot

all_boxplot <- (knz_boxplot | sbc_boxplot) /
               (sev_boxplot | plot_spacer())
all_boxplot


# Looking at differences across datasets ----------------------------------



m1 <- glmmTMB(bray ~ type*dataset + (1|site_year),
              data = all_bray,
              beta_family())

summary(m1)

em <- emmeans(m1, pairwise ~ type | dataset)

em

