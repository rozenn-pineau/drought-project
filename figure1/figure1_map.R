#Makes map for Figure 1

library(gsheet)
library(ggplot2)
library(dplyr)
library(extrafont)
library(colourvalues)
library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
library(MetBrewer)
library(maps)
library(mapdata)
#library(maptools)  #for shapefiles - was removed from CRAN - do I need it?
library(scales)  #for transparency
library("rnaturalearth")
#install.packages("rnaturalearthhires")
library("rnaturalearthhires")
library(rnaturalearth)
library(dplyr)
library(raster)
library(sf)
library(tidyverse)
library(ggrepel)
library(data.table)
library(ggthemes)
library(geodata)
library(ggnewscale)
library(raster)

rm(list= ls())

# ------------------------------------------------------------------------------------- #
# prep data
# ------------------------------------------------------------------------------------- #
#Contemporary samples
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/4.phenotypes")
ld50 <- fread("ld50.txt") 
ld50$long_cols <- color_values(ld50$long, palette = "purples")
ag_idx <- which(ld50$env == "Ag")
nat_idx <- which(ld50$env == "Nat")

# Jitter lat lon for Ag samples so that they are easier to see on the map
ld50$long_mod <- ld50$long
ld50$long_mod[ag_idx] <- ld50$long_mod[ag_idx] - 0.5
ld50$lat_mod <- ld50$lat
ld50$lat_mod[ag_idx] <- ld50$lat_mod[ag_idx] - 0.3

#Herbarium samples
herb <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/5.herbarium/data/herbarium_samps.txt", sep = "\t", header = T)
herb_fil <- herb[herb$Year>=1880,]

herb_fil$cols <- colour_values(herb_fil$Year, palette = "viridis")
min(herb_fil$Long, na.rm=T)



states    <- c("New York","New Jersey","Delaware","Maryland",
               "Pennsylvania","West Virginia","Virginia",
               "North Carolina","South Carolina","Georgia","Florida",
               "Alabama","Mississippi","Louisiana","Arkansas","Texas",
               "Tennessee","Kentucky","Ohio","Indiana","Illinois",
               "Michigan","Wisconsin","Minnesota","Iowa","Missouri",
               "Kansas","Nebraska","Oklahoma",
               "Lake Michigan","Lake Ontario","Lake Superior")
provinces <- c("Ontario")

us <- gadm(country = "USA", level = 1, path = tempdir())
# us_sf <- st_as_sf(us)
# centroids <- st_centroid(us_sf)
#coords <- st_coordinates(centroids)
canada <- gadm(country="CAN",level=1)


us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

#lakes
lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

#ocean
ocean <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'ocean', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

# rivers
rivers <- rnaturalearth::ne_download(scale = 10, 
                                     type = 'rivers_lake_centerlines', 
                                     category = 'physical')  %>% 
  sf::st_as_sf(lakes110, crs = 4269)


#get US states outlines
in_sf <- ne_states(geounit = "United States of America",
                   returnclass = "sf")

uss <- st_as_sf(us.states) %>% 
  mutate(
    lon = map_dbl(geometry, ~st_centroid(.x)[[1]]),
    lat = map_dbl(geometry, ~st_centroid(.x)[[2]]))

#move Michigan label and add Ontario label
uss[uss$NAME_1=="Michigan",]$lon<--85.0554
uss[uss$NAME_1=="Michigan",]$lat<-44.00902
uss<-uss %>% add_row(NAME_1 = "Ontario", lon = -78.5554, lat=45)


# Fortify the GADM objects
#us.states.df <- fortify(us.states)
#ca.provinces.df <- fortify(as.data.frame(ca.provinces))

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/4.phenotypes/")
pdf(file=paste("map_figure1.pdf", sep="" ), bg = "transparent", width=10, height=5)

plain1<- 
  ggplot()+
  geom_sf(data = uss,
          fill = "white",
          color = "black") +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          , fill="#cee5ed", alpha = 0.5)  +
  geom_sf(data = ocean,
          mapping = aes(geometry = geometry),
          fill="#cee5ed", alpha = 0.5)  +
  geom_sf(data = rivers,
          mapping = aes(geometry = geometry),
          color = "#99ccff", alpha=.75) +
  coord_sf(xlim=c(min(herb_fil$Long, na.rm=T),max(herb_fil$Long, na.rm=T)), 
           ylim=c(min(herb_fil$Lat, na.rm=T), 46)) +
  theme_minimal()

#plot map base with contemporary samples
p1 <- plain1 + 
  geom_jitter(data=ld50,
              aes(y = lat_mod, x = long_mod, fill = long_cols, shape = env), size=4, alpha=1, color = "grey50") +
  scale_shape_manual(values = c(Ag = 21, Nat = 22)) +
  scale_fill_manual(values = setNames(ld50$long_cols, ld50$long_cols)) +
  xlab("Longitude") + ylab("Latitude") + 
  theme(legend.position = "none") + 
  new_scale_color() +
  geom_jitter(data=herb_fil, #plot map base with herbarium samples
              aes(y = Lat, x = Long, color = cols), size=2,  alpha=1) + 
  scale_color_manual(values = setNames(herb_fil$cols, herb_fil$cols)) 
p1

dev.off()

# Export legend
pdf(file=paste("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/4.phenotypes/scale_herb_orange.pdf", sep="" ), bg = "transparent", width=6, height=3)

par(family = "Times New Roman", cex = 1.3)
scale_data <- data.frame(a = 4, x = min(herb_fil$Year):max(herb_fil$Year))
scale_data$cols <- colour_values(scale_data$x, palette = "viridis")
barplot(height = scale_data$a, col = scale_data$cols, space = 0, border = NA, names.arg = round(scale_data$x), yaxt='n', ann=FALSE) # 

dev.off()





