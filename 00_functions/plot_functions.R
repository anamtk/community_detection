# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Rhat graph --------------------------------------------------------------


rhat_graph_fun <- function(list, parms, rhat){

  #this creates a dtaaframe out of the Rhat values from the model
  df <- data.frame(id = names(list),
                   Rhat = unlist(lapply(list, paste, collapse = ","))) %>%
    #splits Rhat by , when that list element had more than one value
    mutate(Rhat = str_split(Rhat, ",")) %>%
    #unnests - so makes each Rhat a new row in the df
    unnest(c(Rhat)) %>%
    #make sure Rhat is a numeric
    mutate(Rhat = as.numeric(Rhat)) %>%
    filter(id %in% parms)
  
  #plot histogram and make sure all below 1.1
  plot <- ggplot(df, aes(x = Rhat)) +
    geom_histogram() +
    geom_vline(xintercept = rhat, linetype = 2) +
    theme_bw() +
    scale_y_sqrt() +
    #this is facetted by the parameter ID so you
    # can find problematic parameters
    facet_wrap(~ id)
  
  return(plot)
}


# SAM covariate effects plot function -------------------------------------


effects_plot_fun <- function(model){
  b0 <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(parm == "b0") %>%
    dplyr::select(`50%`) %>%
    as_vector()
  
  betas <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, "b")) %>%
    filter(!str_detect(parm, "b0")) %>%
    filter(!str_detect(parm, 'sig.web'))
  
  beta_plot <- ggplot(betas, aes(x = `50%`, y= parm)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.4) +
    geom_point() +
    labs(x = "Covariate effect", y = "") +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8),
          plot.title.position = "plot",
          panel.grid = element_blank(),
          axis.ticks.y = element_blank())
  
  return(beta_plot)
  
}


# SAM partial plot funciton -----------------------------------------------



partial_plot_fun <- function(model, covariate, df, ID, yearID, start, end, weight, diss){
  
  beta <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(parm == covariate) %>%
    dplyr::select(`50%`) %>%
    as_vector()
  
  b0 <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parm") %>%
    filter(str_detect(parm, "b0")) %>%
    dplyr::select(`50%`) %>%
    summarise(b0 = mean(`50%`, na.rm = T)) %>%
    as_vector()
  
  temp <- df %>%
    dplyr::select(ID, yearID, start:end) %>% #adjust if needed
    pivot_longer(start:end,
                 names_to = "lag",
                 values_to = "var") %>%
    mutate(var = scale(var)) %>%
    pivot_wider(names_from = "lag",
                values_from = "var") %>%
    dplyr::select(-ID, -yearID) %>%
    as.matrix()
  
  #make scaled data long format to get mean and sd
  scale <- df %>%
    dplyr::select(ID, yearID, start:end) %>% #adjust if needed
    pivot_longer(start:end,
                 names_to = "lag",
                 values_to = "var") 
  
  #get mean and SD of OG data to back-transform stuff
  Mean <- mean(scale$var, na.rm = T)
  SD <- sd(scale$var, na.rm = T)
  
  #get weights per lag
  wt <- as.data.frame(model$quantiles) %>%
    rownames_to_column(var = "parameter") %>%
    filter(str_detect(parameter, weight)) %>%
    dplyr::select(`50%`) %>%
    as_vector()
  
  #get tmax dataset
  regT <- df %>%
    dplyr::select(ID, yearID, diss, start:end)
  
  #multiply months by their weights
  regT$Ant <- apply(temp, MARGIN = 1, FUN = function(x){sum(x*wt)})
  
  #revert Tmax to OG data scale
  regT <- regT %>%
    dplyr::select(Ant, diss) %>%
    mutate(Var = Ant*SD + Mean)
  
  #regression prediction for Temperature
  regT <- regT %>%
    mutate(reg = b0 + beta*Ant,
           plogisreg = plogis(reg))
  
  plot <- ggplot(regT) +
    geom_point(aes(x = Var, y = .data[[diss]]), alpha = 0.5,
               position = position_jitter()) +
    geom_line(aes(x = Var, y = plogisreg), size = 1) 
  
  return(plot)
  
}  



