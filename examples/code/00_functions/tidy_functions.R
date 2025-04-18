# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse") 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}



# Extract Chlrophyl for each point ----------------------------------------


extract_chl <- function(files){
  
  # import the list of dataframes in the vector "files", 
  #now a list of dataframes, into R
  #then make them all rasters
  data_list <-  lapply(files, function(x) raster(x))
  
  # make a shortened name of the files
  files_short <- files %>%
    str_replace(".*/", "") %>% # remove the file path
    str_replace("AQUA_MODIS.", "") %>%
    str_replace(".L3m.MO.CHL.chlor_a.4km.nc", "") # remove the file extension
  
  # make the list of files have the name of their file
  names(data_list) <- files_short
  
  data_list2 <- lapply(data_list, function(x) crop(x, ROI))
  
  data_list3 <- lapply(data_list2, function(x) raster::extract(x, extract.pts, sp = T))
 
  return(data_list3)
}


# Make Chl_a values dataframes --------------------------------------------

chla_df <- function(vector){
  
  names <- c("chla", "site")
  df <- as.data.frame(vector) %>%
    dplyr::mutate(site = c("ABUR", "AHND", "AQUE", "BULL", "CARP", "GOLB",
                           "IVEE", "MOHK", "NAPL", "SCDI", "SCTW")) 
  
  colnames(df) <- names
  return(df)
  
}


# Get initials for omega --------------------------------------------------

get_omega_init <- function(cov_matrix){
t_cov <- cov(t1)
t_cov2 <- cov(t2)
t_cov3 <- cov(t3)
#get mean value of diagonal values that are not 0
diag_mean <- mean(diag(t_cov)[diag(t_cov) != 0])

#set all zero values on diagonal to be that mean value
#set all diagonals to the mean (this did work):
diag(t_cov) <- diag_mean


#top and bottom 5% of off-diagonal
(upper <- quantile(t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) > 0)], 
                   probs = c(0.95)))
(lower <- quantile(t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) < 0)], 
                   probs = c(0.05)))

#find mean of values that are positivie but less than the upper quantile
umean <- mean(t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) < upper) & ((t_cov) > 0)])
#set all extreme positive values to this mean
t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) >= upper)] <- umean

#find the mean of values tha are negative but greater than the lower quantile
lmean <- mean(t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) > lower) & ((t_cov) < 0)])
#set all extreme negative values to this mean
t_cov[!(col(t_cov) == row(t_cov)) & ((t_cov) < 0) & (t_cov <= lower)] <- lmean

#get off diagonal mean that is not 0
odiag_mean <- mean(t_cov[!t_cov == 0])
#set any zero values to that off diagonal mean
t_cov[t_cov == 0] <- odiag_mean

#invert to get precision matrix
omega.init <- ginv(t_cov)

}


# Scaled-unscaled dataframe for predictive graphing -----------------------


#get a dataframe to predict onto that includes
# the original scale of the variable 
# the scaled variable, and the "level"
# of the scaled variable. You can
# specify different "lengths", but 20 is a good
# length to choose from 

#x =the original variable
#length = the lenght of the scale (how many values to predict onto)
# name = the name of the variable that you want to attach to the final df
scale_df <- function(x, length, name){
  
  #get the sequenceo f the variable
  var1 <- seq(min(x, na.rm = T),
              max(x, na.rm = T),
              #set the length you want
              length.out = length)
  
  #scale that variable
  varS <- (var1 - mean(x, na.rm = T))/sd(x, na.rm = T)
  
  #make the OG and scaled a DF with appropriate 
  # names and a "level" for graphing later
  final_df <- var1 %>%
    cbind(varS) %>%
    as.data.frame() %>%
    dplyr::rename(!!name := ".") %>%
    dplyr::mutate(level = row_number())
  
  #get the final DF back from the function
  return(final_df)
}



# Rhat DF -----------------------------------------------------------------


rhat_df_fun <- function(list){
  
  #this creates a dtaaframe out of the Rhat values from the model
  df <- data.frame(id = names(list),
                   Rhat = unlist(lapply(list, paste, collapse = ","))) %>%
    #splits Rhat by , when that list element had more than one value
    mutate(Rhat = str_split(Rhat, ",")) %>%
    #unnests - so makes each Rhat a new row in the df
    unnest(c(Rhat)) %>%
    #make sure Rhat is a numeric
    mutate(Rhat = as.numeric(Rhat)) 
  
  return(df)
  
}
