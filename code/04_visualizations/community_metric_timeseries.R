#Exploring yearly community metrics
#Ana Miller-ter Kuile
#July 5, 2023

#this is a script that explores visuals of the enviromentla
#covariates in different ways...

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 

package.list <- c("here", "tidyverse",
                  "data.table", "patchwork",
                  "viridis")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

## And loading them

for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())
# Load data ---------------------------------------------------------------

data <- read.csv(here("data_outputs",
                         "community_stability",
                         "stability_metrics_with_covariates.csv"))

# Explore -----------------------------------------------------------------

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data$gain_d <- get_density(data$YEAR, data$gain, n = 100)
data$loss_d <- get_density(data$YEAR, data$loss, n = 100)
data$temp_d <- get_density(data$YEAR, data$TEMP_C, n = 100)


data %>%
  summarise(gain = mean(gain),
            loss = mean(loss),
            Temp = mean(TEMP_C, na.rm = T),
            kelp = mean(DRY_GM2))

(a <- ggplot(data, aes(x = YEAR, y = gain, color = gain_d)) +
  geom_hline(yintercept = 0.153466, linetype = 2) +
  geom_point() +
  scale_color_viridis() +
  scale_x_continuous(breaks = c(2001:2023)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_vline(xintercept = 2013.5, linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = 2016.5, linetype = 2, alpha = 0.4)  )


b <- ggplot(data, aes(x = YEAR, y = loss, color = loss_d)) +
  geom_hline(yintercept = 0.1548491, linetype = 2) +
  geom_point() +
  scale_color_viridis() +
  scale_x_continuous(breaks = c(2001:2023)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_vline(xintercept = 2013.5, linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = 2016.5, linetype = 2, alpha = 0.4)

c <- ggplot(data, aes(x = YEAR, y = TEMP_C)) +
  geom_hline(yintercept = 16.30706, linetype = 2) +
  geom_point() +
  scale_x_continuous(breaks = c(2001:2023)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_vline(xintercept = 2013.5, linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = 2016.5, linetype = 2, alpha = 0.4)

(timeseries1 <- a / c / b)

d <- ggplot(data, aes(x = YEAR, y = chla)) +
  geom_hline(yintercept = 16.30706, linetype = 2) +
  geom_point() +
  scale_x_continuous(breaks = c(2001:2023)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_vline(xintercept = 2013.5, linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = 2016.5, linetype = 2, alpha = 0.4)

(timeseries2 <- a / d / b)


ggplot(data, aes(x = YEAR, y = DRY_GM2)) +
  geom_hline(yintercept = 324.0753, linetype = 2) +
  geom_point() +
  scale_x_continuous(breaks = c(2001:2023)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_vline(xintercept = 2013.5, linetype = 2, alpha = 0.4) +
  geom_vline(xintercept = 2016.5, linetype = 2, alpha = 0.4)
