
# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "jagsUI",
                  'rjags',
                  'mcmcplots',
                  "coda",
                  'patchwork') #mcmc output


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

data_list <- readRDS(here('02_konza_birds',
                          "data_outputs",
                          'SAM',
                          "model_inputs",
                          "knz_bray_SAM_input_data_raw.RDS"))


# Parameters to save ------------------------------------------------------

params <- c('b',
            'b0',
            'wA',
            'wB',
            #'wC',
            'var.process')



# JAGS model --------------------------------------------------------------

model <- here('02_konza_birds',
              "code", 
              "02_kb_analyses",
              'SAM',
              'raw',
              "jags",
              "bird_SAM_raw.R")

Sys.time()
mod <- jagsUI::jags(data = data_list,
                    inits = NULL,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 3,
                    n.iter = 2835,
                    n.burnin = 1500,
                    DIC = TRUE)

Sys.time()

# Check convergence -------------------------------------------------------

mcmcplot(mod$samples)

gelman.diag(mod$samples, multivariate = F)

rhats <- mod$Rhat

parm <- c("b", "b0", "wA", "wB", "var.process", "deviance")

rhat_graph_fun(list =rhats, parms = parm, rhat = 1.1) +
  labs(title = "KNZ LTER bird SAM Rhat")

ggsave(plot = last_plot(),
       filename = here("pictures",
                       "supplementary",
                       'SAM',
                       "KNZ_SAM_raw_Rhat_graph.jpg"),
       height = 4,
       width = 6,
       units = "in")
# Output summary ----------------------------------------------------------

sum <- summary(mod$samples)

saveRDS(sum, here('02_konza_birds',
                  'data_outputs',
                  'SAM',
                  'model_outputs',
                  'knz_SAM_summary_raw.RDS'))


# Update for GOF ----------------------------------------------------------
params2 <- c('bray.rep',
            'resid')

modGOF <- update(mod,
                 parameters.to.save = params2,
                 n.iter = 1335)

sumGOF <- summary(modGOF$samples)

saveRDS(sumGOF, here('02_konza_birds',
                  'data_outputs',
                  'SAM',
                  'model_outputs',
                  'knz_SAM_raw_GOF_summary.RDS'))

# Check interaction for overfitting ---------------------------------------


# Interaction -------------------------------------------------------------

# #this interaction is overfitting the data, so i removed it from the model

# bird_bray <- read.csv(here('02_konza_birds',
#                            'data_outputs',
#                            "SAM",
#                            'data_prep',
#                            'stability_metrics_with_covariates.csv'))
# 

# blTP <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[4]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# blT <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[1]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# blP <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(parm == "b[2]") %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# b0 <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parm") %>%
#   filter(str_detect(parm, "b0")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# #get temparutres on scaled scale
# temp_temp <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, TAVE:TAVE_l5) %>% #adjust if needed
#   pivot_longer(TAVE:TAVE_l5,
#                names_to = "lag",
#                values_to = "temp") %>%
#   mutate(temp = scale(temp)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "temp") %>%
#   dplyr::select(-WATERSHED, -RECYEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# tmaxscale <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, TAVE:TAVE_l5) %>% #adjust if needed
#   pivot_longer(TAVE:TAVE_l5,
#                names_to = "lag",
#                values_to = "temp")
# 
# #get mean and SD of OG data to back-transform stuff
# mean <- mean(tmaxscale$temp, na.rm = T)
# sd <- sd(tmaxscale$temp, na.rm = T)
# 
# #get weights per month
# t_wt <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wA")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# 
# #get kelp on scaled scale
# ppt_temp <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, PPT:PPT_l5) %>% #adjust if needed
#   pivot_longer(PPT:PPT_l5,
#                names_to = "lag",
#                values_to = "ppt") %>%
#   mutate(ppt = scale(ppt)) %>%
#   pivot_wider(names_from = "lag",
#               values_from = "ppt") %>%
#   dplyr::select(-WATERSHED, -RECYEAR) %>%
#   as.matrix()
# 
# #make scaled data long format to get mean and sd
# pptscale <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, PPT:PPT_l5) %>% #adjust if needed
#   pivot_longer(PPT:PPT_l5,
#                names_to = "lag",
#                values_to = "ppt") 
# 
# #get mean and SD of OG data to back-transform stuff
# meanp <- mean(pptscale$ppt, na.rm = T)
# sdp <- sd(pptscale$ppt, na.rm = T)
# 
# #get weights per month
# p_wt <- as.data.frame(sum$quantiles) %>%
#   rownames_to_column(var = "parameter") %>%
#   filter(str_detect(parameter, "wB")) %>%
#   dplyr::select(`50%`) %>%
#   as_vector()
# 
# #get tmax dataset
# regT <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, bray, TAVE:TAVE_l5)
# 
# #multiply months by their weights
# regT$TAnt <- apply(temp_temp, MARGIN = 1, FUN = function(x){sum(x*t_wt)})
# 
# #revert Tmax to OG data scale
# regT <- regT %>%
#   dplyr::select(TAnt, bray, WATERSHED, RECYEAR) %>%
#   mutate(Temp = TAnt*sd + mean)
# 
# #kelp dataset
# regP <- bird_bray %>%
#   dplyr::select(WATERSHED, RECYEAR, bray, PPT:PPT_l5)
# 
# #multiply months by their weights
# regP$PAnt <- apply(ppt_temp, MARGIN = 1, FUN = function(x){sum(x*p_wt)})
# 
# #revert Tmax to OG data scale
# regP <- regP %>%
#   dplyr::select(PAnt, bray, WATERSHED, RECYEAR) %>%
#   mutate(PPT = PAnt*sd + mean)
# 
# #regression
# regB <- regT %>%
#   left_join(regP, by = c("WATERSHED", "RECYEAR", "bray"))
# 
# temp2 <- scale_df(x = regB$Temp,
#                   length = 20,
#                   name = "temp")
# 
# ppt2 <- scale_df(x = regB$PPT,
#                   length = 20,
#                   name = "ppt") %>%
#   rename("varP" = "varS",
#          "levelP" = "level")
# 
# tp <- temp2 %>%
#   cross_join(ppt2) %>%
#   mutate(reg = b0 + blT*varS + blP*varP + blTP*varS*varP,
#          plogisreg = plogis(reg))
# 
# (a <- ggplot(tp, aes(x = temp, y = ppt, fill = plogisreg)) +
#   geom_tile() +
#   geom_contour(aes(z = plogisreg), color = "white") +
#   scale_fill_viridis_c() +
#   theme(axis.title = element_blank()))
# 
# (b <- ggplot(regB, aes(x = PPT)) +
#   geom_boxplot() +
#   coord_flip() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()))
# 
# (c <- ggplot(regB, aes(x = Temp)) +
#   geom_boxplot() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()))
# 
# (b + a)/(plot_spacer() + c) +
#   plot_layout(widths = c(1, 3),
#               heights = c(3, 1))
# 



