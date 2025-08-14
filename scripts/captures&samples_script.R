# Created: August 6, 2025
# Updated: August 14, 2025
# 
# Author: Carina Isabella Motta
# 
# Captures & Samples

# 1. LOAD PACKAGES -------------------------------------------------------------

# a vector listing package names needed 

package.list <- c("here", #so I don't have to deal with setting a WD
                  "vegan", #species accumulation curves
                  "tidyverse", #data cleaning
                  "dplyr" #data cleaning
)

#creating another list of new packages (if there are any)
new.packages <- package.list[!(package.list %in% installed.packages()
                               [,"Package"])]

#installing the packages if they aren't already on the computer
if(length(new.packages)) install.packages(new.packages)

#and loading the packages into R with a for loop
for(i in package.list){library(i, character.only = T)}

# 2. LOAD DATA  ----------------------------------------------------------------
## 2.1 Capture data----

# captures
captures <- readr::read_csv(here::here("data", "bird_captures.csv"))

# corresponding family names 
bird.spp.fam <- readr::read_csv(here::here("data", "bird_spp_family.csv"))

## 2.2 Droppings and seed sample data------

#seed count and identification 
seeds <- readr::read_csv(here::here("data", "droppings_seeds.csv"))

#seed morphospecies 
morpho <- readr::read_csv(here::here("data", "droppings_seeds_morphotypes.csv"))

# 3. CLEAN AND PREP DATA ####
## 3.1 Capture data ------

#transform into a tibble to facilitate using packages like dyplr 
captures <- as_tibble(captures)

#cut excess rows
captures <- captures[-(492:999),]

#cut excess columns
captures <- captures[,-(28:36)]

#subset captures to only include certain columns: capture date, session, day, 
#plot, time the net was opened and time net was closed, species, sample number, 
#whether feces were collected, if the bird is a recap  (y = 1/n = 0), and if 
#the bird is fruit-eating

subset.caps <- captures %>% 
  select(1:6, 11:12, 19, 23, 25, 26)

#cut excess column
bird.spp.fam <- bird.spp.fam[,-(3)]

#add to the capture dataset
subset.caps <- subset.caps %>%
  left_join(bird.spp.fam, by = "bird_species") %>%
  relocate(bird_family, .after = nets_closed) 

## 3.2 Sample data----

#join seed count and morphospecies data 

#transform into a tibble to facilitate using packages like dyplr 
morpho <- as_tibble(morpho)

#cut excess rows
morpho <- morpho[-(42:987),]

#select columns of interest
subset.morpho <- morpho %>% 
  select(1:5)

#select columns of interest
subset.seeds <- seeds %>% 
  select(1:7)

#merge seed and morphospp data 
seed.morpho <- merge(x=subset.seeds, y=subset.morpho,
                     by="seed_morpho", all.x =T)

## 3.3 Merge capture and sample data-----

#append seed data to subset.caps
subset.caps$seeds <- ifelse(subset.caps$sample %in% seed.morpho$sample, 1, 0)


seeds.sample <- seed.morpho %>%
  group_by(sample) %>%
  summarise(no_seeds = sum(no_seeds),
            no_morpho_spp = max(no_morpho))

captures.seeds <-  subset.caps %>%
  left_join(seeds.sample, by = "sample")

#fill NAs with 0s
captures.seeds$no_seeds[is.na(captures.seeds$no_seeds)] <- 0

#fill NAs with 0s
captures.seeds$no_morpho_spp[is.na(captures.seeds$no_morpho_spp)] <- 0

#reorder the columns
captures.seeds <- captures.seeds[,c(1:9, 11:16, 10)]

#write.csv(captures.seeds, here::here("data", "processed_data",
#                                     "captures.seeds.csv"))


# 4. INDIVIDUAL CAPTURE SUMMARY & METRICS ####

## 4.1 Captures, feces samples, & seeds collected------

#filter out the no capture days
caps.summary <- captures.seeds[!grepl("no_captures", 
                                      captures.seeds$bird_species),]

caps.total.summary <- caps.summary %>%
  summarise(n_caps = n_distinct(sample),
            n_recaps = sum(recapture == 1, na.rm = T),
            n_caps_frug = sum(frugivorous == 1, na.rm = T),
            n_spp_birds = n_distinct(bird_species),
            n_spp_frug = n_distinct(bird_species[frugivorous == 1]),
            n_samples = sum(feces_collected == "y", na.rm = T),
            n_samples_frug = sum(feces_collected == "y" & frugivorous == 1),
            n_samples_seeds = sum(seeds == 1, na.rm = T),
            n_seeds = sum(no_seeds)
            
  )

## 4.2 Net hours-----

net.hours <- captures.seeds %>%
  distinct(date, .keep_all = TRUE) %>%
  mutate(
    nets_opened = as.POSIXct(nets_opened, format = "%H:%M:%S"), 
    nets_closed = as.POSIXct(nets_closed, format = "%H:%M:%S"),
    net_hours = as.numeric(difftime(nets_closed, nets_opened, 
                                    units = "hours")) * 5
  )

# sum total net hours
net.hours.total <- sum(net.hours$net_hours)

# Print result
print(net.hours.total) #2763

caps.total.summary$net_hours <- sum(net.hours$net_hours)

caps.total.summary$cap_rate <- caps.total.summary$n_caps/caps.total.summary$net_hours

caps.total.summary$cap_rate_frug <- caps.total.summary$n_caps_frug/caps.total.summary$net_hours

print(caps.total.summary, width = Inf)

# 5. SPECIES-LEVEL SUMMARY & METRICS ####

## 5.1 Captures, feces samples, & seeds collected----- 

#create a list of species captured, with the number of captures, number of samples, and the mass of the species 
caps.spp.summary <- caps.summary %>%
  group_by(bird_species, bird_family, frugivorous) %>%
  summarise(
    no_caps = n(),  # Counts the number of captures for each species
    no_samples = sum(feces_collected == "y", na.rm = T), 
    mass_caps = mean(mass, na.rm = TRUE)
  ) %>%
  ungroup()

#summarize samples with seeds and number of seeds per species to create a species-level summary table 
seeds.summary <- subset.seeds %>%
  group_by(bird_species) %>%
  summarise(
    no_samples_seeds = n_distinct(sample),  
    no_seed_total = sum(no_seeds) 
  ) %>%
  ungroup()

#merge seeds summary and species capture summary 
caps.spp.seed.summary <- merge(x=caps.spp.summary, y=seeds.summary,
                               by="bird_species", all.x =T)

#fill NAs with 0s
caps.spp.seed.summary[is.na(caps.spp.seed.summary)] <- 0

#reorder the columns
species.summary <- caps.spp.seed.summary[,c(2, 1, 4, 5, 7, 8, 3, 6)]

## 5.2 Capture rate-----

#calculate the capture rate by dividing the number of captures by the total number of net hours 
species.summary$cap_rate <- (species.summary$no_caps/2763) 

write.csv(species.summary, here::here("data", "processed_data",
                                       "species.summary.csv"))

# 6. PLOT-LEVEL SUMMARY & METRICS ####

## 6.1 Captures, feces samples, & seeds collected----

#number of captures, recaps per plot, number of frugivorous captures, number of samples, number of samples from frugivores, number of samples with seeds
caps.plot.summary <- caps.summary %>% #use cap.summary, that already filtered out no capture days
  group_by(plot) %>%
  summarise(no_captures = n_distinct(sample),
            no_recaps = sum(recapture == 1, na.rm = T),
            no_caps_frug = sum(frugivorous == 1, na.rm = T),
            no_spp_birds = n_distinct(bird_species),
            no_spp_frug = n_distinct(bird_species[frugivorous == 1]),
            no_samples = sum(feces_collected == "y", na.rm = T),
            no_samples_frug = sum(feces_collected == "y" & frugivorous == 1),
            no_samples_seeds = sum(seeds == 1, na.rm = T)
            
  )


## 6.2 Net hours----

# calculate net hours per plot
net.hours.plot <- net.hours %>%
  group_by(plot) %>%  # Group by plot
  summarise(total_net_hours = round(sum(net_hours)))  # Sum net hours per plot  

# append to plot.summary 
plot.summary <- merge(x=caps.plot.summary, y=net.hours.plot,
                      by="plot", all.x =T) 

## 6.3 Capture rate-----

#all captures
plot.summary$cap_rate <- (plot.summary$no_caps/plot.summary$total_net_hours)

#frugivorous captures
plot.summary$cap_rate_frug <- (plot.summary$no_caps_frug/plot.summary$total_net_hours)

cor(plot.summary[, c("cap_rate", "cap_rate_frug")])

write.csv(plot.summary, here::here("data", "processed_data",
                                    "plot.summary.csv"))

