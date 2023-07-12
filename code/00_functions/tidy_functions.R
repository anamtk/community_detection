

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
