#SAM model results
#Ana Miller-ter Kuile
#October 27, 2023

#this script generates results figures for SAM models

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

source(here('00_functions',
            'tidy_functions.R'))
# Load data ---------------------------------------------------------------

fish_sam <- readRDS(here('01_sbc_fish',
                         'monsoon',
                         'SAM',
                         'outputs',
                         'fish_SAM_summary.RDS'))

fish_bray <- read.csv(here("01_sbc_fish",
                          "data_outputs",
                          'SAM',
                          'data_prep',
                          "stability_metrics_with_covariates.csv"))

sev_sam <- readRDS(here('03_sev_grasshoppers',
                        'monsoon',
                        'SAM',
                        'outputs',
                        'sev_SAM_summary.RDS'))

sev_bray <- read.csv(here("03_sev_grasshoppers",
                          "data_outputs",
                          "SAM",
                          'data_prep',
                          'sev_stability_metrics_with_covariates.csv'))


# Effect plots ------------------------------------------------------------

#fish
b0 <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b0") %>%
  dplyr::select(`50%`) %>%
  as_vector()

bs <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "b")) %>%
  filter(!str_detect(parm, "b0"))

beta_sbc <- ggplot(bs, aes(x = parm, y = `50%`)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  labs(x = "Covariate", y = "Covariate effect \n (Median and 95% BCI)",
       title = "SBC fish") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  scale_x_discrete(labels = c("Kelp Biomass", "Temperature")) +
  coord_flip() +
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 15),
        plot.title.position = "panel",
        plot.title = element_text(hjust = 0.5))

#grasshoppers
b0sev <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b0") %>%
  dplyr::select(`50%`) %>%
  as_vector()

bsev <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "b")) %>%
  filter(!str_detect(parm, "b0")) %>%
  filter(parm != "sig.web")

beta_sev <- ggplot(bsev, aes(x = parm, y = `50%`)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  labs(x = "Covariate", y = "Covariate effect \n (Median and 95% BCI)",
       title = "SEV grasshoppers") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
  scale_x_discrete(labels = c("Temperature", "Plant biomass")) +
  coord_flip()+
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 15),
        plot.title.position = "panel",
        plot.title = element_text(hjust = 0.5))

beta_sev / beta_sbc

ggsave(plot = last_plot(),
       filename = here('pictures',
                       'sam_models',
                       'sam_covariate_effects.jpg'),
       height = 7, 
       width = 6,
       units = "in")
# Partial plots -----------------------------------------------------------


# Fish partial plots ------------------------------------------------------

### temperature
blT <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b[2]") %>%
  dplyr::select(`50%`) %>%
  as_vector()

b0 <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "b0")) %>%
  dplyr::select(`50%`) %>%
  summarise(b0 = mean(`50%`, na.rm = T)) %>%
  as_vector()

#get temparutres on scaled scale
temp_temp <- fish_bray %>%
  dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
  pivot_longer(TEMP_C:TEMP_C_l5,
               names_to = "lag",
               values_to = "temp") %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = "lag",
              values_from = "temp") %>%
  dplyr::select(-SITE_TRANS, -YEAR) %>%
  as.matrix()

#make scaled data long format to get mean and sd
tmaxscale <- fish_bray %>%
  dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
  pivot_longer(TEMP_C:TEMP_C_l5,
               names_to = "lag",
               values_to = "temp") 

#get mean and SD of OG data to back-transform stuff
mean <- mean(tmaxscale$temp, na.rm = T)
sd <- sd(tmaxscale$temp, na.rm = T)

#get weights per month
t_wt <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(str_detect(parameter, "wB")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get tmax dataset
regT <- fish_bray %>%
  dplyr::select(SITE_TRANS, YEAR, bray, TEMP_C:TEMP_C_l5)

#multiply months by their weights
regT$TAnt <- apply(temp_temp, MARGIN = 1, FUN = function(x){sum(x*t_wt)})

#revert Tmax to OG data scale
regT <- regT %>%
  dplyr::select(TAnt, bray) %>%
  mutate(Temp = TAnt*sd + mean)

#regression prediction for Temperature
regT <- regT %>%
  mutate(reg = b0 + blT*TAnt,
         plogisreg = plogis(reg))

(fisht <- ggplot(regT) +
  geom_point(aes(x = Temp, y = bray), alpha = 0.5) +
  geom_line(aes(x = Temp, y = plogisreg), size = 1) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity",
       title = "SBC fish") +
  #annotate(geom = "text", x = 13.25, y = 0.4, label = "More similar") + 
  #annotate(geom = "text", x = 13.25, y = 0.8, label = "More different") +
  theme(plot.title.position = "panel",
        plot.title = element_text(hjust = 0.5)))
  
###WEIGHTS
tweights <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]") ~ "Warm",
                            parm %in% c("wB[2]", "wB[4]", "wB[6]") ~ "Cold")) %>%
  mutate(year = case_when(parm == "wB[1]" ~ 0,
                          parm %in% c("wB[2]", 'wB[3]') ~ 1,
                          parm %in% c('wB[4]', 'wB[5]') ~ 2,
                          parm == 'wB[6]' ~ 3))

warmcol <- '#d8b365'
coldcol <- '#5ab4ac'

(fish_tweights <- tweights %>%
  ggplot(aes(x = year, y= `50%`, shape = season)) +
  geom_hline(yintercept = 1/6, linetype = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                position = position_dodge(width = 0.5),
                width = 0) +
  scale_color_manual(values = c(Warm = warmcol, Cold = coldcol)) +
  labs(x = "Years into the past",
       y = "Importance weight\n(median and 95% BCI)"))

fishtgraphs <- fisht + fish_tweights
# Interaction -------------------------------------------------------------

#this interaction is overfitting the data, so i removed it from the model
# blTK <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[3]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector() 
# 
# ### temperature
# blT <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[2]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# b0 <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b0")) %>%
#   dplyr::select(`50%`) %>%
#   summarise(b0 = mean(`50%`, na.rm = T)) %>%
#   as_vector()
# 
# #get temparutres on scaled scale
# temp_temp <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
#   pivot_longer(TEMP_C:TEMP_C_l5,
#                names_to = "lag",
#                values_to = "temp") %>%
#   mutate(temp = scale(temp)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "temp") %>%
#   dplyr::select(-SITE_TRANS, -YEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# tmaxscale <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, TEMP_C:TEMP_C_l5) %>% #adjust if needed
#   pivot_longer(TEMP_C:TEMP_C_l5,
#                names_to = "lag",
#                values_to = "temp") 
# 
# #get mean and SD of OG data to back-transform stuff
# mean <- mean(tmaxscale$temp, na.rm = T)
# sd <- sd(tmaxscale$temp, na.rm = T)
# 
# #get weights per month
# t_wt <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wB")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# 
# #get kelp on scaled scale
# kelp_temp <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>% #adjust if needed
#   pivot_longer(DRY_GM2:DRY_GM2_l5,
#                names_to = "lag",
#                values_to = "kelp") %>%
#   mutate(kelp = scale(kelp)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "kelp") %>%
#   dplyr::select(-SITE_TRANS, -YEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# kelpscale <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, DRY_GM2:DRY_GM2_l5) %>% #adjust if needed
#   pivot_longer(DRY_GM2:DRY_GM2_l5,
#                names_to = "lag",
#                values_to = "kelp") 
# 
# #get mean and SD of OG data to back-transform stuff
# meank <- mean(kelpscale$kelp, na.rm = T)
# sdk <- sd(kelpscale$kelp, na.rm = T)
# 
# #get weights per month
# k_wt <- as.data.frame(fish_sam$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wA")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# #get tmax dataset
# regT <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, bray, TEMP_C:TEMP_C_l5)
# 
# #multiply months by their weights
# regT$TAnt <- apply(temp_temp, MARGIN = 1, FUN = function(x){sum(x*t_wt)})
# 
# #revert Tmax to OG data scale
# regT <- regT %>%
#   dplyr::select(TAnt, bray, SITE_TRANS, YEAR) %>%
#   mutate(Temp = TAnt*sd + mean)
# 
# #kelp dataset
# regK <- fish_bray %>%
#   dplyr::select(SITE_TRANS, YEAR, bray, DRY_GM2:DRY_GM2_l5)
# 
# #multiply months by their weights
# regK$KAnt <- apply(kelp_temp, MARGIN = 1, FUN = function(x){sum(x*k_wt)})
# 
# #revert Tmax to OG data scale
# regK <- regK %>%
#   dplyr::select(KAnt, bray, SITE_TRANS, YEAR) %>%
#   mutate(Kelp = KAnt*sd + mean)
# 
# #regression
# regB <- regT %>%
#   left_join(regK, by = c("SITE_TRANS", "YEAR", "bray")) 
# 
# temp2 <- scale_df(x = regB$Temp,
#                   length = 20,
#                   name = "temp")
#   
# kelp2 <- scale_df(x = regB$Kelp,
#                   length = 20,
#                   name = "kelp") %>%
#   rename("varK" = "varS",
#          "levelK" = "level")
# 
# tk <- temp2 %>%
#   cross_join(kelp2) %>%
#   mutate(reg = b0 + blT*varS + blTK*varS*varK,
#          plogisreg = plogis(reg))
# 
# a <- ggplot(tk, aes(x = temp, y = kelp, fill = plogisreg)) +
#   geom_tile() +
#   geom_contour(aes(z = plogisreg), color = "white") +
#   scale_fill_viridis_c() +
#   theme(axis.title = element_blank())
# 
# b <- ggplot(regB, aes(x = Kelp)) +
#   geom_boxplot() +
#   coord_flip() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# 
# c <- ggplot(regB, aes(x = Temp)) +
#   geom_boxplot() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# 
# (b + a)/(plot_spacer() + c) +
#   plot_layout(widths = c(1, 3),
#               heights = c(3, 1))


# Grasshopper partial plots -----------------------------------------------

#temperature
### temperature
blTs <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b[1]") %>%
  dplyr::select(`50%`) %>%
  as_vector()

b0s <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "b0")) %>%
  dplyr::select(`50%`) %>%
  summarise(b0 = mean(`50%`, na.rm = T)) %>%
  as_vector()

#get temparutres on scaled scale
temp_temps <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, Temp:Temp_l5) %>% #adjust if needed
  pivot_longer(Temp:Temp_l5,
               names_to = "lag",
               values_to = "temp") %>%
  mutate(temp = scale(temp)) %>%
  pivot_wider(names_from = "lag",
              values_from = "temp") %>%
  dplyr::select(-site_web_trans, -YEAR) %>%
  as.matrix()

#make scaled data long format to get mean and sd
tscalesev <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, Temp:Temp_l5) %>% #adjust if needed
  pivot_longer(Temp:Temp_l5,
               names_to = "lag",
               values_to = "temp") 

#get mean and SD of OG data to back-transform stuff
means <- mean(tscalesev$temp, na.rm = T)
sds <- sd(tscalesev$temp, na.rm = T)

#get weights per month
t_wts <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(str_detect(parameter, "wA")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get tmax dataset
regTs <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, bray, Temp:Temp_l5)

#multiply months by their weights
regTs$TAnt <- apply(temp_temps, MARGIN = 1, FUN = function(x){sum(x*t_wts)})

#revert Tmax to OG data scale
regTs <- regTs %>%
  dplyr::select(TAnt, bray) %>%
  mutate(Temp = TAnt*sds + means)

#regression prediction for Temperature
regTs <- regTs %>%
  mutate(reg = b0s + blTs*TAnt,
         plogisreg = plogis(reg))

(hoppert <- ggplot(regTs) +
  geom_point(aes(x = Temp, y = bray), 
             position = position_jitter(width = 0.2),
             alpha = 0.5) +
  geom_line(aes(x = Temp, y = plogisreg), size = 1) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity",
       title = "SEV grasshoppers") +
  annotate(geom = "text", x = 12.75, y = 0.3, label = "More similar") + 
  annotate(geom = "text", x = 12.75, y = 0.7, label = "More different")+
    theme(plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5)))

###WEIGHTS
tweightss <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wA")) %>%
  mutate(season = case_when(parm %in% c("wA[1]", "wA[3]", "wA[5]") ~ "Cold",
                            parm %in% c("wA[2]", "wA[4]", "wA[6]") ~ "Warm")) %>%
  mutate(year = case_when(parm %in% c("wA[1]", "wA[2]") ~ 0,
                          parm %in% c('wA[3]', 'wA[4]') ~ 1,
                          parm %in% c('wA[5]', 'wA[6]') ~ 2))

warmcol <- '#d8b365'
coldcol <- '#5ab4ac'

(hopper_tweights <- tweightss %>%
    ggplot(aes(x = year, y= `50%`, shape = season)) +
    geom_hline(yintercept = 1/6, linetype = 2) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                  position = position_dodge(width = 0.5),
                  width = 0) +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    scale_color_manual(values = c(Warm = warmcol, Cold = coldcol)) +
    labs(x = "Years into the past",
         y = "Importance weight\n(median and 95% BCI)"))

hoppertgraphs <- hoppert + hopper_tweights

#NPP
blNs <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b[2]") %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get temparutres on scaled scale
temp_npp <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, NPP:NPP_l10) %>% #adjust if needed
  pivot_longer(NPP:NPP_l10,
               names_to = "lag",
               values_to = "npp") %>%
  mutate(npp = scale(npp)) %>%
  pivot_wider(names_from = "lag",
              values_from = "npp") %>%
  dplyr::select(-site_web_trans, -YEAR) %>%
  as.matrix()

#make scaled data long format to get mean and sd
nscale <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, NPP:NPP_l10) %>% #adjust if needed
  pivot_longer(NPP:NPP_l10,
               names_to = "lag",
               values_to = "npp") 

#get mean and SD of OG data to back-transform stuff
meannpp <- mean(nscale$npp, na.rm = T)
sdnpp <- sd(nscale$npp, na.rm = T)

#get weights per month
n_wts <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parameter") %>%
  filter(str_detect(parameter, "wC")) %>%
  dplyr::select(`50%`) %>%
  as_vector()

#get tmax dataset
regNs <- sev_bray %>%
  dplyr::select(site_web_trans, YEAR, bray, NPP:NPP_l10)

#multiply months by their weights
regNs$NAnt <- apply(temp_npp, MARGIN = 1, FUN = function(x){sum(x*n_wts)})

#revert Tmax to OG data scale
regNs <- regNs %>%
  dplyr::select(NAnt, bray) %>%
  mutate(NPP = NAnt*sdnpp + meannpp)

#regression prediction for Temperature
regNs <- regNs %>%
  mutate(reg = b0s + blNs*NAnt,
         plogisreg = plogis(reg))

(hoppern <- ggplot(regNs) +
    geom_point(aes(x = NPP, y = bray), 
               position = position_jitter(width = 0.2),
               alpha = 0.5) +
    geom_line(aes(x = NPP, y = plogisreg), size = 1) +
    labs(x = "Plant biomass",
         y = "Bray-Curtis Dissimilarity") )
    #annotate(geom = "text", x = 11.25, y = 0.3, label = "More similar") + 
    #annotate(geom = "text", x = 11.25, y = 0.7, label = "More different"))

###WEIGHTS
nweights <- as.data.frame(sev_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wC")) %>%
  mutate(season = case_when(parm %in% c("wC[1]", "wC[3]", "wC[5]", "wC[7]",
                                        'wC[9]', 'wC[11]') ~ "Fall",
                            parm %in% c("wC[2]", "wC[4]", "wC[6]",
                                        'wC[8]', 'wC[10]') ~ "Spring")) %>%
  mutate(year = case_when(parm %in% c("wC[1]", "wC[2]") ~ 0,
                          parm %in% c('wC[3]', 'wC[4]') ~ 1,
                          parm %in% c('wC[5]', 'wC[6]') ~ 2,
                          parm %in% c('wC[7]', 'wC[8]') ~ 3,
                          parm %in% c('wC[9]', 'wC[10]') ~ 4,
                          parm %in% c('wC[11]') ~ 5))

warmcol <- '#d8b365'
coldcol <- '#5ab4ac'

(hopper_nweights <- nweights %>%
    ggplot(aes(x = year, y= `50%`, shape = season)) +
    geom_hline(yintercept = 1/12, linetype = 2) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                  position = position_dodge(width = 0.5),
                  width = 0) +
    scale_x_continuous(breaks = c(0, 1,2,3,4,5)) +
    scale_color_manual(values = c(Warm = warmcol, Cold = coldcol)) +
    labs(x = "Years into the past",
         y = "Importance weight\n(median and 95% BCI)"))

hoppernppgraphs <- hoppern + hopper_nweights


hoppertgraphs / hoppernppgraphs / fishtgraphs


ggsave(filename = here('pictures',
                       'sam_models',
                       'sam_partial_plots.jpg'),
       height = 7, 
       width = 8,
       units = "in")




