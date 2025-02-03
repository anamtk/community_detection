#functional diversity calculations
#Ana Miller-ter Kuile
#January 27, 2025

#this script calculates functional diversity of 
#each iteration of the bird dataset MSAM

#Fric, Rao

#https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html

# Load packages -----------------------------------------------------------

package.list <- c("here", "tidyverse",
                  'readxl', 'patchwork',
                  'FD', 'fundiversity', 'ape')

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())

# Data --------------------------------------------------------------------

birds <- read.csv(here('examples',
                       "data_output",
                       'bird_fundiv',
                       'tidy_data',
                       'bird_msam_tidy_data.csv')) %>%
  mutate(AOUCODE = case_when(COMMONNAME == "American Goldfinch" ~ 'AGOL',
                             TRUE ~ AOUCODE))

codes <- read.csv(here('examples',
                       'data_raw',
                       'bird_fundiv',
                       'IBP-AOS-LIST23.csv'))

avonet3 <- read.csv(here('examples',
                         'data_raw',
                         'bird_fundiv',
                         'AVONET_eBird.csv'))

# Prep trait data ---------------------------------------------------------

#bird size
konz_codes <- unique(birds$AOUCODE)

#get the codes from IBP that only 
#match Konza birds
traits <- codes %>%
  filter(SPEC %in% konz_codes) %>%
  dplyr::select(-SP, -CONF, -SPEC6, -CONF6) %>%
  left_join(avonet3, by = c("SCINAME" = "Species2")) %>%
  dplyr::select(SPEC, COMMONNAME, SCINAME,
                #diet traits:
                Trophic.Level,
                Trophic.Niche,
                #behavioral:
                Primary.Lifestyle,
                Migration,
                Habitat, 
                Habitat.Density,
                #morphology
                Mass, 
                #Tarsus.Length,
                Tail.Length,
                Wing.Length,
                #Kipps.Distance, 
                Hand.Wing.Index,
                Beak.Width, 
                #Beak.Depth,
                #Beak.Length_Nares, 
                Beak.Length_Culmen
  )



# Prep link between model species ID and species  -------------------------

bird_ids <- birds %>%
  distinct(COMMONNAME, AOUCODE, SpecID)

traits2 <- traits %>%
  left_join(bird_ids, by = c("COMMONNAME", "SPEC" = "AOUCODE")) %>%
  dplyr::select(COMMONNAME, Trophic.Level:Beak.Length_Culmen) 

#

# species x traits matrix -------------------------------------------------

#need a species x traits matrix
#all continuous traits need to be standardized
#with categorical and continuous traits, need to do
#two more steps:
## compute gower distance
## run PCOA and pull out trait values from the first
## axis that explains variation as "traits" data

# Get gower dist of triats matrix -----------------------------------------

#continuous variables need to be scaled
#categorical need to be factored and then numerical
trait_matdf <- traits2 %>%
  mutate(Habitat.Density = as.factor(Habitat.Density),
         Migration = as.factor(Migration),
         Trophic.Level = as.factor(Trophic.Level)) %>%
  mutate(Habitat = factor(Habitat, 
                          levels = c("Human Modified", "Forest",
                                     "Grassland","Shrubland",     
                                     "Rock", "Wetland", 
                                     "Riverine","Woodland" )),
         Trophic.Niche = factor(Trophic.Niche,
                                levels = c("Granivore","Nectarivore",  
                                           "Frugivore", "Omnivore",
                                           "Invertivore",'Vertivore',
                                           "Aquatic predator")),
         Primary.Lifestyle = as.factor(Primary.Lifestyle)) %>%
  mutate_if(is.factor, as.numeric)#

trait_mat <- trait_matdf %>%
  column_to_rownames(var = 'COMMONNAME') %>%
  mutate_if(is.numeric, scale) %>%
  as.matrix()

trait_gower <- FD::gowdis(trait_mat)

pcoa <- ape::pcoa(trait_gower)

#use these as the matrix for the speciesxtraits
spp_pcoa_scores <- pcoa$vectors %>% 
  as_tibble(rownames = NA) %>%  
  rownames_to_column("Species") 

#selecting first five for now as a test
pcoa_mat <- pcoa$vectors[,1:5] 


# Site x species matrix ---------------------------------------------------


#get the site x species matrix
site_sp_matrix <- birds %>%
  left_join(bird_ids, by = c("SpecID", "COMMONNAME", "AOUCODE")) %>%
  dplyr::select(COMMONNAME, TransID, yrID, NOBS) %>%
  group_by(COMMONNAME, TransID, yrID) %>%
  summarise(maxN = max(NOBS, na.rm = T)) %>%
  ungroup() %>%
  unite("trans_yr",
        c("TransID", "yrID"),
        sep = "_") %>%
  pivot_wider(names_from = "COMMONNAME",
              values_from = "maxN") %>%
  column_to_rownames(var = 'trans_yr') %>%
  as.matrix()


# Try... ------------------------------------------------------------------

#accounts for abundances
rao <- fd_raoq(pcoa_mat, 
               site_sp_matrix)

rao2 <- rao %>%
  separate(site,
           into = c("TransID", "yrID"))

#just presence-absence
fric <- fundiversity::fd_fric(pcoa_mat, 
                              site_sp_matrix, 
                              stand = T)

fric2 <- fric %>%
  separate(site,
           into = c("TransID", "yrID"))

FD <- rao2 %>%
  left_join(fric2, by = c("TransID", "yrID"))

saveRDS(FD, here('examples',
                 "data_output",
                 'bird_fundiv',
                 'tidy_data',
                 'bird_fd_metrics_raw.RDS'))


# Modeled -----------------------------------------------------------------

samples <- readRDS(here('examples',
                        "data_output",
                        'bird_fundiv',
                        'computing_cluster_outputs',
                        'bird_N_samples.RDS'))

# Prep modeled N dataset --------------------------------------------------

ndf <- as.data.frame.table(samples$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         SpecID = as.numeric(as.factor(Var2)),
         TransID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, SpecID, TransID, yrID, N)

# Site x species matrix ---------------------------------------------------

fd_fun <- function(iter){
  
  #get the site x species matrix
  site_sp_matrix <- ndf %>%
    filter(Iteration == {{iter}}) %>%
    left_join(bird_ids, by = c("SpecID")) %>%
    dplyr::select(-AOUCODE, -SpecID, -Iteration) %>%
    filter(!is.na(N)) %>%
    unite("trans_yr",
          c("TransID", "yrID"),
          sep = "_") %>%
    pivot_wider(names_from = "COMMONNAME",
                values_from = "N") %>%
    column_to_rownames(var = 'trans_yr') %>%
    as.matrix()
  
  #accounts for abundances
  rao <- fd_raoq(pcoa_mat, 
                 site_sp_matrix) %>%
    separate(site,
             into = c("TransID", "yrID")) 
  
  #just presence-absence
  fric <- fundiversity::fd_fric(pcoa_mat, 
                                site_sp_matrix, 
                                stand = T)%>%
    separate(site,
             into = c("TransID", "yrID"))
  
  FD <- rao %>%
    left_join(fric, by = c("TransID", "yrID")) %>%
    mutate(Iteration = {{iter}})
  
  return(FD)
  
}

iter_vec <- ndf %>%
  distinct(Iteration) %>%
  as_vector()

iter_vec1 <- iter_vec[1:50]

fd_list <- lapply(iter_vec1, fd_fun)

fd_df <- bind_rows(fd_list)

iter_vec2 <- iter_vec[51:100]

fd_list2 <- lapply(iter_vec2, fd_fun)

fd_df2 <- bind_rows(fd_list2)

iter_vec3 <- iter_vec[101:150]

fd_list3 <- lapply(iter_vec3, fd_fun)

fd_df3 <- bind_rows(fd_list3)

iter_vec4 <- iter_vec[151:200]

fd_list4 <- lapply(iter_vec4, fd_fun)

fd_df4 <- bind_rows(fd_list4)

iter_vec5 <- iter_vec[201:250]

fd_list5 <- lapply(iter_vec5, fd_fun)

fd_df5 <- bind_rows(fd_list5)

iter_vec6 <- iter_vec[251:300]

fd_list6 <- lapply(iter_vec6, fd_fun)

fd_df6 <- bind_rows(fd_list6)

iter_vec7 <- iter_vec[301:350]

fd_list7 <- lapply(iter_vec7, fd_fun)

fd_df7 <- bind_rows(fd_list7)

iter_vec8 <- iter_vec[351:400]

fd_list8 <- lapply(iter_vec8, fd_fun)

fd_df8 <- bind_rows(fd_list8)

iter_vec9 <- iter_vec[401:451]

fd_list9 <- lapply(iter_vec9, fd_fun)

fd_df9 <- bind_rows(fd_list9)

iter_vec10 <- iter_vec[452:500]

fd_list10 <- lapply(iter_vec10, fd_fun)

fd_df10 <- bind_rows(fd_list10)

fd_all <- bind_rows(fd_df, fd_df2, fd_df3, fd_df4,
                    fd_df5, fd_df6, fd_df7,
                    fd_df8, fd_df9, fd_df10)

# Summarise FD ------------------------------------------------------------

FD_sum <- fd_all %>%
  group_by(TransID, yrID) %>%
  summarise(Q_mean = mean(Q, na.rm = T),
            Q_sd = sd(Q, na.rm = T),
            FRic_mean = mean(FRic, na.rm = T),
            FRic_sd = sd(FRic, na.rm = T))

saveRDS(FD_sum, here('examples',
                     "data_output",
                     'bird_fundiv',
                     'tidy_data',
                     'bird_fd_metrics_corrected.RDS'))

