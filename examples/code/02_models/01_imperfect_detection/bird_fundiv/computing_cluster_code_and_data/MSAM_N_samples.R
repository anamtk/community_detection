#Get N out of converged model
#Ana Miller-ter Kuile
#September 19, 2023

#this script pulls out samples for N and calculates 
#functional diversity from them

# Load packages ---------------------------------------------------------------
Sys.time()


# Load packages,
package.list <- c('jagsUI',"coda",'dplyr', 
                  'stringr','magrittr', 'purrr',
                  'tidyr','ggplot2', 'forcats', 'tidyr',
                  'tibble','FD', 'fundiversity', 'ape') 


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}


# Load model --------------------------------------------------------------

mod <- readRDS( file ="/scratch/atm234/konza_birds/outputs/bird_MSAM_RE_model2.RDS")

# Update model to track bray ----------------------------------------------

parms <- c("N")

mod2 <- update(mod,
               parameters.to.save = parms,
               parallel = TRUE,
               n.iter = 5000)

# Get the N samples -------------------------------------------------------

samples <- mod2$sims.list

saveRDS(samples, 
        file = "/scratch/atm234/konza_birds/outputs/bird_N_samples.RDS")

# Load other datasets -----------------------------------------------------

birds <- read.csv(file = "/scratch/atm234/konza_birds/inputs/knz_tidy_data_for_model.csv")

birds <- birds %>%
  mutate(AOUCODE = case_when(COMMONNAME == "American Goldfinch" ~ 'AGOL',
                             TRUE ~ AOUCODE))

codes <- read.csv(file = "/scratch/atm234/konza_birds/inputs/IBP-AOS-LIST23.csv")

avonet3 <- read.csv(file = "/scratch/atm234/konza_birds/inputs/AVONET_eBird.csv")


# Prep modeled N dataset --------------------------------------------------

ndf <- as.data.frame.table(samples$N) %>%
  mutate(Iteration = as.numeric(as.factor(Var1)),
         SpecID = as.numeric(as.factor(Var2)),
         TransID = as.numeric(as.factor(Var3)),
         yrID = as.numeric(as.factor(Var4))) %>%
  rename(N = Freq) %>%
  dplyr::select(Iteration, SpecID, TransID, yrID, N)

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
site_sp_matrix <- ndf %>%
  left_join(bird_ids, by = c("SpecID")) %>%
  dplyr::select(-AOUCODE, -SpecID) %>%
  filter(!is.na(N)) %>%
  unite("it_trans_yr",
        c("Iteration", "TransID", "yrID"),
        sep = "_") %>%
  pivot_wider(names_from = "COMMONNAME",
              values_from = "N") %>%
  column_to_rownames(var = 'it_trans_yr') %>%
  as.matrix()


# Try... ------------------------------------------------------------------

#accounts for abundances
rao <- fd_raoq(pcoa_mat, 
               site_sp_matrix)

saveRDS(rao, file = "/scratch/atm234/konza_birds/outputs/bird_rao.RDS")

#just presence-absence
fric <- fundiversity::fd_fric(pcoa_mat, 
                              site_sp_matrix, 
                              stand = T)

saveRDS(fric, file = "/scratch/atm234/konza_birds/outputs/bird_fric.RDS")

