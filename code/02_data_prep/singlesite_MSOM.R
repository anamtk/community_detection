# Ana Miller-ter Kuile
# January 17, 2023



# Load packages -----------------------------------------------------------


package.list <- c("here", "tidyverse")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load data ---------------------------------------------------------------

benthic <- read.csv(here("data_raw", 
                         "Annual_Cover_All_Years_20220809.csv"))




# Clean benthic dataset ---------------------------------------------------


colnames(benthic)

#2 species codes for
#Stephanocystis osmundacea - CH and CYOS - this is holdfast vs. young
#Macrocystis pyrifera - DMH and MH - this is holdfast vs. young
#Amphiroa beauvoisii - AMZO, UEC - this is two spp code with no explanation...

#remove unidentified species from the group
benthic1 <- benthic %>%
  filter(!str_detect(SCIENTIFIC_NAME, "Unidentified|Unidentifiable")) %>%
  #one species code is confusingly "NA" so have to tell computer
  # that isn't missing data
  mutate(SP_CODE = case_when(is.na(SP_CODE) ~ "NA",
                             TRUE ~ SP_CODE)) %>%
  group_by(SCIENTIFIC_NAME, SITE, TRANSECT, QUAD, SIDE, YEAR) %>%
  #some species have more than one code, so combining here (might
  # change this appraoch later)
  summarise(PERCENT_COVER1 = sum(PERCENT_COVER)) %>%
  ungroup()

benthic2 <- benthic1 %>%
  #get occupancy from percent cover
  mutate(OCC = case_when(PERCENT_COVER1 > 0 ~ 1,
                         PERCENT_COVER1 == 0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  #get unique side_quad combos
  unite(SIDE_QUAD, c("QUAD", "SIDE"), sep = "_") %>%
  #get unique site_transect combos
  unite(SITE_TRANSECT, c("SITE", "TRANSECT"), sep = "_") %>%
  #select variables of interest
  dplyr::select(YEAR,SIDE_QUAD, SITE_TRANSECT, SCIENTIFIC_NAME, OCC) %>%
  #remove missing occupancy values now 
  #might want to set to zero? Think about this later down the road
  filter(!is.na(OCC)) 


# Subset a Transect for single site model ---------------------------------

#subset one of the transects for the single-site model
benthic_test <- benthic2 %>%
  filter(SITE_TRANSECT == "AQUE_1") 

write.csv(benthic_test, 
          here("data_outputs",
               "raw_community",
               "single_site_MSOM_long.csv"))

# #a function that makes the site, year, occupancy data
# # into a matrix of species by year
# #filled with occupacy for that species for that year
# df_function <- function(df, quad){
#   matrix <- df %>%
#     #select the qudrat of choice
#     filter(SIDE_QUAD == quad) %>%
#     #make years columns and cells filled with occupancy
#     #for a psecies for that year
#     pivot_wider(names_from = "YEAR",
#                 values_from = "OCC") %>%
#     #get rid of the variables no longer needed
#     dplyr::select(-SITE_TRANSECT, -SIDE_QUAD) %>%
#     #make the species names the column names
#     column_to_rownames('SCIENTIFIC_NAME') %>%
#     mutate(across(.cols = everything(), as.double)) %>%
#     #make it a matrix
#     as.matrix()
#   
#   return(matrix) #return the matrix from the function
# }
# 
# #a list of the four quadrats in that transect so we can get
# # our replicates
# list <- c("20_I", "20_O", "40_I", "40_O")
# 
# #using the list above, make a list of dataframes
# #usig the DF function on the benthic_test df and 
# #mapping/looping trhough the list (quad) argument
# ytemp <- map(list, 
#                   ~ df_function(df = benthic_test,
#                                 quad = .))
# 
# y <- array(unlist(ytemp), dim = c(122,22,4))

y <- benthic_test %>%
  group_by(YEAR, SCIENTIFIC_NAME) %>%
  summarise(occ = sum(OCC, na.rm = T)) %>%
  pivot_wider(names_from = "SCIENTIFIC_NAME",
              values_from = "occ") %>%
  mutate(across(everything(), ~replace_na(.,0))) %>%
  mutate(across(everything(), as.integer)) %>%
  column_to_rownames(var = "YEAR") %>%
  as.matrix()

write.csv(y, 
          here("data_outputs",
               "raw_community",
               "single_site_MSOM_matrix.csv"))

z <- (y>0)*1
z[z == 0] <- NA
# Get covariates and list elements ----------------------------------------

n.species <- ncol(y)
n.years <- nrow(y)
n.reps <- 4

visibility <- read.csv(here("data_raw",
                            "Annual_fish_comb_20220809.csv"))

visibility1 <- visibility %>%
  filter(SITE == "AQUE" & TRANSECT == 1) %>%
  distinct(YEAR, DATE, SITE, TRANSECT, VIS)  %>%
  filter(VIS != -99999)

vis <- visibility1 %>%
  dplyr::select(VIS) %>%
  mutate(VIS = scale(VIS)) %>%
  as_vector()


# Export all data objects for model ---------------------------------------

data <- list(y = y,
             n.species = n.species,
             n.years = n.years,
             n.reps = n.reps,
             vis = vis,
             z = z)

saveRDS(data, file = here("data_outputs",
                          "model_inputs",
                          "data_singlesite_MSOM.RDS"))
