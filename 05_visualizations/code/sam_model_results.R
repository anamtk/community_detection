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


# Partial plots -----------------------------------------------------------

b0 <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(parm == "b0") %>%
  dplyr::select(`50%`) %>%
  as_vector()

bs <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "b")) %>%
  filter(!str_detect(parm, "b0"))

ggplot(bs, aes(x = parm, y = `50%`)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  scale_x_discrete(labels = c("Kelp", "Temperature", "Kelp*Temp")) +
  coord_flip()


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

fisht <- ggplot(regT) +
  geom_point(aes(x = Temp, y = bray), alpha = 0.5) +
  geom_line(aes(x = Temp, y = plogisreg), linewidth = 1) +
  labs(x = "Temperature",
       y = "Bray-Curtis Dissimilarity") +
  annotate(geom = "text", x = 13.25, y = 0.4, label = "More similar") + 
  annotate(geom = "text", x = 13.25, y = 0.8, label = "More different")
  
###WEIGHTS
tweights <- as.data.frame(fish_sam$quantiles) %>%
  rownames_to_column(var = "parm") %>%
  filter(str_detect(parm, "wB")) 

warmcol <- '#d8b365'
coldcol <- '#5ab4ac'

fish_tweights <- tweights %>%
  mutate(season = case_when(parm %in% c("wB[1]", "wB[3]", "wB[5]") ~ "Warm",
                            parm %in% c("wB[2]", "wB[4]", "wB[6]") ~ "Cold")) %>%
ggplot(aes(x = parm, y= `50%`, shape = season)) +
  geom_hline(yintercept = 1/6, linetype = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
  scale_x_discrete(labels = c("0", "1",
                              "2", "3", "4",
                              "5")) +
  scale_color_manual(values = c(Warm = warmcol, Cold = coldcol)) +
  labs(x = "Seasons into the past",
       y = "Importance weight\n(median and 95% BCI)")

fisht + fish_tweights
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










