## Created: 02-October-2024
## Updated: 06-June-2025

## Author: Carina Isabella Motta

## Evaluating a  DNA metabarcoding approach of bird feces to quantify and 
## predict seed dispersal

# Data cleaning and analysis to compare two data sets: 
# Morphological identifications of seeds and 
# identification of plant DNA in same bird feces 

# 1 LOAD PACKAGES------------------------------------------------------------

# a vector listing package names needed 

package.list <- c("here", #makes this project the working directory
                  "tidyverse", #data cleaning
                  "dplyr", #data cleaning
                  "reshape2", #matrix manipulation
                  "VennDiagram", #VennDiagram figure
                  "ggalluvial" #interaction web 
)


#installing the packages if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()
                               [,"Package"])]

if(length(new.packages)) install.packages(new.packages)

#and loading the packages into R with a for loop
for(i in package.list){library(i, character.only = T)}

# 2 LOAD DATA---------------------------------------------------------------

## 2.1 LOAD SEED DATA----

seeds <- readr::read_csv(here::here("data", 
                                    "droppings_seeds.csv"))

## 2.2 LOAD SEED MORPHOSPECIES DATA----

morpho <- readr::read_csv(here::here("data",
                                     "droppings_seeds_morphotypes.csv"))

## 2.3 LOAD METABARCODING DATA----

gen <- readr::read_csv(here::here("data",
                                  "droppings_metabarcoding.csv"))


# 3 DATA CLEANING AND PREP-----------------------------------------------------

## 3.1 SEED & SEED MORPHOLOGY DATA----

### 3.1.1 SEED DATA-----

# subset columns that we are interested in 
seeds <- seeds %>% 
  select(3, 4, 6, 7, 8)

### 3.1.2 MORPHOLOGY DATA-----

# subset columns that we are interested in 
morpho <- morpho %>% 
  select(1:5)

### 3.1.3 MERGE SEED AND SEED MORPHOLOGY DATA----

seeds.morph <- merge(x=seeds, y=morpho, by="seed_morpho", all.x =T)

##3.2 METABARCODING DATA----

# cut excess rows and columns 

#cut excess rows
gen <- gen[-(561:1246),]

#cut excess columns 
gen <- gen[,-(16:32)]

#filter to only include pools that include samples with seeds, and corresponding blanks

gen.seeds <- gen %>%
    filter(pool %in% c("P_00", "P_01", "P_02", "P_03","P_04", "P_05", "P_06",
                       "P_07", "P_08", "P_09","P_10", "P_11", "P_12"))%>%
    filter(seeds == 1 | grepl("blank", sample))


# 4 SEED & SEED MORPHOLOGY SUMMARY STATISTICS----
## 4.1 ALL SAMPLES----

n_distinct(seeds.morph$unique)

n_distinct(seeds.morph$plant_family)

n_distinct(seeds.morph$plant_genus)

sum(seeds.morph$no_seeds)

bird.spp.all <- seeds.morph %>%
                  group_by(bird_species) %>%
                  summarize(n_all = n_distinct(sample))

## 4.2 SUCESSFULLY SEQUENCED SAMPLES----
gen.seeds.seq <- gen.seeds %>% 
                 filter(final_fate == 1)

#cross reference seed morphology data set with sequencing to only include 
#the samples that were successfully sequenced. That way we can create our 
#interaction webs using only data that is present in both data sets

seeds.morph.seq <- seeds.morph %>%
  filter(sample %in% gen.seeds.seq$sample)

#check to see if was correctly filtered, should have 52 samples now 

n_distinct(seeds.morph.seq$sample) # 52, correctly filtered 

#check the number of morphospecies present after filtering out the samples 
#that didn't successfully sequence, should be 36

n_distinct(seeds.morph.seq$unique) # 36 morphospecies

n_distinct(seeds.morph.seq$plant_family) #21 families

# calculate number of genera, excluding morphospecies identified to family or
# unknown

n_distinct(seeds.morph.seq %>%
    filter(!plant_genus == "NA") %>%
    pull(plant_genus)
)

# number of seeds 
sum(seeds.morph.seq$no_seeds)

### 4.2.1 MAKE TABLE 1----

#create a table of how many samples of each species was successfully 
#sequenced 

bird.spp.successful <- seeds.morph.seq %>%
  group_by(bird_species) %>%
  summarize(n_successful = n_distinct(sample))
# create a table of how many samples of each species were collected vs.
# successfully sequenced (Table 1)

bird.spp.samples <- merge(x=bird.spp.all, y=bird.spp.successful, 
                          by="bird_species", all.x =T) 

#fill NAs with 0s
bird.spp.samples[is.na(bird.spp.samples)] <- 0

# save table as a .csv
write.csv(bird.spp.samples, here::here("Analyses",
                                       "figures&results", 
                                       "Table_1.csv"))


### 4.2.2 CALCULATE NUMBER OF PLANT TAXA PER SAMPLE ----

seeds.summary.ind <- seeds.morph.seq %>% 
  group_by(sample) %>%
  summarise(n_taxa = n_distinct(plant_genus))

mean(seeds.summary.ind$n_taxa)

sd(seeds.summary.ind$n_taxa)

### 4.2.3 CALCULATE NUMBER OF SAMPLES PER PLANT TAXA ----
seeds.summary.plants <- seeds.morph.seq %>% 
  group_by(plant_genus) %>%
  summarise(n_samples = n_distinct(sample))

mean(seeds.summary.plants$n_samples)

sd(seeds.summary.plants$n_samples)

## 4.3 SEED MORPHOLOGY TAXA SUMMARY (APPENDIX S1)-----------------------------

seed.summary <- seeds.morph.seq %>%
  group_by(plant_genus_morpho) %>%
  summarise(
    plant_family = first(plant_family),
    plant_genus = first(plant_genus),
    n_droppings = n_distinct(sample), 
    n_bird_spp = n_distinct(bird_species), 
    seed_sum = sum(no_seeds),
    .groups = "drop"
  )

write.csv(seed.summary, here::here("analyses",
                                   "figures&results", 
                                   "Appendix_S1.csv"))

# 5 METABARCODING SUMMARY STATISTICS----
## 5.1 FATE SUMMARIES----

# final fates of each OTU
# 0 means it didn't make the read threshold,
# 1 means is was successful and kept
# 2 means it didn't make the match threshold
# 3 means it wasn't successfully sequenced 
# 4 means it was excluded due to suspected contamination 
# 5 means that it was a blank sample that generated an OTU

sample.fate <- gen.seeds %>%
  filter(!grepl("blank", sample)) %>%
  group_by(sample) %>%
  summarise(
    final_fate = ifelse(any(final_fate == 1), "successful",
                  ifelse(any(final_fate == 0), "insufficient_reads",
                    ifelse(any(final_fate == 2), "insufficient_match",
                      ifelse(any(final_fate == 3),"unsuccessful_seq", 
                        ifelse(any(final_fate == 4), "contamination"))))
    ))

#here we summarize, samples only fell into three categories:
#   successful, unsuccessfully sequenced, or had insufficient reads
#   while individual OTUs were deleted due to insufficient match or 
#   contamination, those two were not motives for excluding samples
sample.fate.summary <- sample.fate %>%
  group_by(final_fate) %>%
  summarise(total_samples = n())

# now we summarize the OTU fates

OTU.fate <- gen.seeds %>%
  filter(!grepl("blank", sample)) %>%
  group_by(OTU) %>%
  summarise(
    final_fate = ifelse(any(final_fate == 1), "successful",
                    ifelse(any(final_fate == 0), "insufficient_reads",
                      ifelse(any(final_fate == 2), "insufficient_match",
                        ifelse(any(final_fate == 3),"unsuccessful_seq", 
                          ifelse(any(final_fate == 4), "contamination", NA))))     ), 
    no_reads_OTU = max(no_reads_OTU))%>%
  mutate(no_reads_OTU = replace_na(no_reads_OTU, 0))

# create summary of OTU fates
OTU.fate.summary <- OTU.fate %>%
  group_by(final_fate) %>%
  summarise(total_OTUs = n_distinct(OTU), n_reads = sum(no_reads_OTU))

## 5.2 READS & OTUs before quality filtering ---- 
t.reads <- gen.seeds %>%
  filter(!grepl("blank", sample),!is.na(OTU),!is.na(no_reads_OTU))  %>%
  group_by(OTU) %>%
  summarise(t_reads = sum(no_reads_OTU))

#the sum of reads
sum(t.reads$t_reads)

#the mean number of reads per OTU
mean(t.reads$t_reads)

#the standard deviation of reads
sd(t.reads$t_reads)

# determine the number of OTUs recovered per sample
t.OTU.sample <- gen.seeds %>%
  filter(!grepl("blank", sample),!is.na(OTU), !is.na(no_reads_OTU))  %>%
  group_by(sample) %>%
  summarise(n_OTU = n_distinct(OTU))

#mean number of OTUs recovered per sample
mean(t.OTU.sample$n_OTU)

#SD of number of OTUs recovered per sample
sd(t.OTU.sample$n_OTU)

## 5.3 READS & OTUs after quality filtering ----- 
t.reads.filtered <- gen.seeds %>%
  filter(final_fate == 1)  %>%
  group_by(OTU) %>%
  summarise(t_reads_filtered = sum(no_reads_OTU))

#the sum of reads
sum(t.reads.filtered$t_reads_filtered)

#the mean number of reads per OTU
mean(t.reads.filtered$t_reads_filtered)

#the standard deviation of reads
sd(t.reads.filtered$t_reads_filtered)

#reads per dropping sample
t.reads.sample.filtered <- gen.seeds %>%
  filter(final_fate == 1)  %>%
  group_by(sample) %>%
  summarise(t_reads_sample_filtered = sum(no_reads_OTU))

mean(t.reads.sample.filtered$t_reads_sample_filtered)

sd(t.reads.sample.filtered$t_reads_sample_filtered)

#OTUs per dropping sample
t.OTU.sample.filtered <- gen.seeds %>%
  filter(final_fate == 1)  %>%
  group_by(sample) %>%
  summarise(n_OTU_filtered = n_distinct(OTU))

mean(t.OTU.sample.filtered$n_OTU_filtered)

sd(t.OTU.sample.filtered$n_OTU_filtered)

## 5.4 PLANT GENERA----

n.plants <- gen.seeds %>%
  filter(!grepl("blank", sample), !is.na(OTU)) %>%
  filter(final_fate == 1)  %>%
  group_by(plant_family, plant_genus) %>%
  summarise(n_plant_fam = n_distinct(sample), 
            n_plant_genera = n_distinct(sample))

n_distinct(n.plants$plant_family)

n_distinct(n.plants$plant_genus)

n_distinct(n.plants %>%
          filter(!plant_genus %in% c("Serjania/Paullinia", "Lippia/Lantana")) %>%
             pull(plant_genus)
)



## 5.5 METABARCODING TAXA SUMMARY (APPENDIX S2)-------------------------------

gen.summary <- gen.seeds.seq %>%
group_by(plant_genus) %>%
  summarise(
    plant_family = first(plant_family),
    n_droppings = n_distinct(sample), 
    n_bird_spp = n_distinct(bird_species), 
    .groups = "drop"
  )

write.csv(gen.summary, here::here("analyses",
                                   "figures&results", 
                                   "Appendix_S2.csv"))
# 6 TRIPARTITE VISUALIZATION ----
## 6.1 CREATE INTERACTION MATRIX WITH GENETIC DATA-----

# Create genetic matrix
gen.matrix <- gen.seeds.seq %>%
  group_by(sample, plant_genus) %>%
  summarise(occurrences = n(), .groups = 'drop') %>%
  pivot_wider(names_from = plant_genus, 
              values_from = occurrences, 
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

# Convert to presence/absence (1/0)
gen.matrix[gen.matrix > 1] <- 1

# Prepare for visualization


# Add rownames as column for melting
gen.matrix.w.rowname <- rownames_to_column(gen.matrix, var = "Sample")

# Transform to long format
genetic.list <- melt(gen.matrix.w.rowname)

#rename columns
colnames(genetic.list) <- c("Sample", "Taxa", "Detections")

# Filter to only positive detections
df.genetic <- subset(genetic.list, Detections >= 1)

# Define factor orders

# Sample order
desired.order.samples <- c(
  "040", "050", "052", "084", "070", "043", "053", "272", 
  "432", "073", "419", "322", "468", "437", "286", "099",
  "110", "136", "287", "257", "412", "381", "113", "447", 
  "321", "103", "246", "269", "253", "132", "271", "384", 
  "181", "241", "242", "243", "064", "212", "255", "482",
  "201", "199", "198", "481", "410", "087", "244", "290", 
  "027", "477", "195", "297"
)

# Genus order
desired.order.taxa <- c(
  "Casearia", "Ocotea", "Terminalia", "Trichilia", "Cestrum",
  "Chamissoa", "Piper", "Cecropia", "Serjania/Paullinia",
  "Centrolobium", "Myriopus", "Tilesia", "Ruprechtia", 
  "Gouania", "Citharexylum", "Solanum", "Syzygium",
  "Banisteriopsis", "Albizia", "Aspidosperma", "Lippia/Lantana",
  "Psychotria", "Alchornea", "Nicotiana", "Amphilophium",
  "Cordia", "Miconia", "Trema", "Senna", "Maclura",
  "Dendropanax", "Stylogyne", "Allophylus", "Eugenia",
  "Guarea", "Syagrus", "Myrciaria", "Myrcia", "Psidium",
  "Morus", "Rubus", "Siparuna"
)

# Set factor levels
df.genetic$Sample <- factor(df.genetic$Sample, levels = desired.order.samples)
df.genetic$Taxa <- factor(df.genetic$Taxa, levels = desired.order.taxa)

# Create alluvial plot

genetic.plot <- ggplot(df.genetic,
                   aes(axis1 = Sample, axis2 = Taxa, y = Detections)) +
  geom_alluvium(aes(fill = Taxa), width = 0.1, alpha = 0.7) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Morus" = "#192609", "Casearia" = "#577882",
      "Trichilia" = "#5C613B", "Miconia" = "#9F4147",
      "Cestrum" = "#AF8932", "Myrcia" = "#DD971D",
      "Nicotiana" = "#D47D11", "Solanum" = "#A2430C",
      "Alchornea" = "#710A07", "Cecropia" = "#6C3A2D",
      "Psychotria" = "#6E7B5A", "Rubus" = "#6EAF7A",
      "Allophylus" = "#5A8C4E", "Maclura" = "#466A21",
      "Psidium" = "#405C1D", "Piper" = "#4B6A50",
      "Myriopus" = "#3A4322", "Tilesia" = "#6C6C7D",
      "Siparuna" = "#855662", "Stylogyne" = "#827B47",
      # All black entries grouped together
      "Ocotea" = "#000000", "Terminalia" = "#000000",
      "Chamissoa" = "#000000", "Serjania/Paullinia" = "#000000",
      "Centrolobium" = "#000000", "Ruprechtia" = "#000000",
      "Gouania" = "#000000", "Citharexylum" = "#000000",
      "Syzygium" = "#000000", "Banisteriopsis" = "#000000",
      "Albizia" = "#000000", "Aspidosperma" = "#000000",
      "Lippia/Lantana" = "#000000", "Amphilophium" = "#000000",
      "Cordia" = "#000000", "Trema" = "#000000",
      "Senna" = "#000000", "Dendropanax" = "#000000",
      "Eugenia" = "#000000", "Guarea" = "#000000",
      "Syagrus" = "#000000", "Myrciaria" = "#000000"
    ))

genetic.plot

# Save outputs

ggsave(
  filename = here::here("analyses", "figures&results", "gen_graph.svg"),
  plot = genetic.plot,
  width = 4,
  height = 16,
  dpi = 720
)

## 6.2 CREATE INTERACTION MATRIX WITH MORPHOLOGY DATA-----

# Create morphology matrix

seed.matrix <- seeds.morph.seq %>%
  mutate(genus_final = ifelse(is.na(plant_genus), plant_genus_morpho, plant_genus)) %>%
  group_by(sample, genus_final) %>%
  summarise(occurrences = n(), .groups = 'drop') %>%
  pivot_wider(names_from = genus_final, 
              values_from = occurrences,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

# Prepare for visualization

# Add rownames as column for melting
seed.matrix.w.rowname <- rownames_to_column(seed.matrix, var = "Sample")

# Transform to long format
seed.list <- melt(seed.matrix.w.rowname)

colnames(seed.list) <- c("Sample", "Taxa", "Detections")

# Ensure presence/absence (1/0)
seed.list$Detections[seed.list$Detections > 1] <- 1


# Filter to only positive detections
df.seed <- subset(seed.list, Detections >= 1)

# Define factor orders

# Genus order for morphology data
desired.order.seed <- c(
  "Casearia", "Trichilia", "unknown_3", "Piper", 
  "unknown_4", "Cestrum", "Chamissoa", "Amaranthaceae",
  "Nicotiana", "Tilesia", "Alchornea", "Cecropia", 
  "unknown_2", "Psychotria", "Myriopus", "Lantanta",
  "Solanum", "Miconia", "Allophylus", "Maclura", 
  "Myrtaceae_2", "Stylogyne", "Fabaceae_1", "Myrcia", 
  "unknown_1", "Psidium", "Morus", "Rubus", "Struthanthus",
  "Pera", "Siparuna"
)

# Set factor levels (using same sample order as genetic data)
df.seed$Sample <- factor(df.seed$Sample, levels = desired.order.samples)
df.seed$Taxa <- factor(df.seed$Taxa, levels = desired.order.seed)


# Create alluvial plot

morph.plot <- ggplot(df.seed,
                     aes(axis1 = Taxa, axis2 = Sample, y = Detections)) +
  geom_alluvium(aes(fill = Taxa), width = 0.1, alpha = 0.7) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Morus" = "#192609",
      "Casearia" = "#577882",
      "Trichilia" = "#5C613B",
      "Miconia" = "#9F4147",
      "Cestrum" = "#AF8932",
      "Myrcia" = "#DD971D",
      "Nicotiana" = "#D47D11",
      "Solanum" = "#A2430C",
      "Alchornea" = "#710A07",
      "Cecropia" = "#6C3A2D",
      "Psychotria" = "#6E7B5A",
      "Rubus" = "#6EAF7A",
      "Allophylus" = "#5A8C4E",
      "Maclura" = "#466A21",
      "Psidium" = "#405C1D",
      "Piper" = "#4B6A50",
      "Myriopus" = "#3A4322",
      "Tilesia" = "#6C6C7D",
      "Siparuna" = "#855662",
      "Stylogyne" = "#827B47",
      # All black entries grouped together
      "Fabaceae_1" = "black",
      "unknown_4" = "black",
      "Pera" = "black",
      "Struthanthus" = "black",
      "Lantanta" = "black",
      "Myrtaceae_2" = "black",
      "unknown_1" = "black",
      "unknown_2" = "black",
      "Amaranthaceae" = "black",
      "unknown_3" = "black"
    )
  )

morph.plot

# save outputs 
ggsave(
  filename = here::here("analyses", "figures&results", "morph_graph.svg"),
  plot = morph.plot,
  width = 4,
  height = 16,
  dpi = 720
)

# 7 JACCARD INDEX ----
## 7.1 CREATE MATRICES OF THE TWO DATASETS WITH MATCHING DIMENSIONS----

# Show ordered summary of genetic detections
gen.list.summary <- genetic.list %>%
  group_by(Taxa) %>%
  summarise(n = sum(Detections))

gen.list.summary.ordered <- gen.list.summary[order(-gen.list.summary$n), ]

print(gen.list.summary.ordered, n = 41)

## 7.2 ALIGN MATRIX DIMENSIONS BY HANDLING MISSING DATA ----

# FIND AND HANDLE DATA MISSING FROM MORPHOLOGY DATASET
missing.from.morph <- genetic.list %>% 
  anti_join(seed.list, by = "Taxa")

missing.from.morph$Detections <- 0

seed.complete <- bind_rows(seed.list, missing.from.morph)

# FIND AND HANDLE DATA MISSING FROM GENETIC DATASET
missing.from.gen <- seed.list %>% 
  anti_join(genetic.list, by = "Taxa")

missing.from.gen$Detections <- 0

gen.complete <- bind_rows(genetic.list, missing.from.gen)

# RECONSTRUCT MATRICES WITH COMPLETE DIMENSIONS
genetic.matrix.complete <- gen.complete %>%
  pivot_wider(names_from = Taxa, values_from = Detections) %>%
  column_to_rownames(var = "Sample")

morpho.matrix.complete <- seed.complete %>%
  pivot_wider(names_from = Taxa, values_from = Detections) %>%
  column_to_rownames(var = "Sample")

# Standardize column order
genetic.matrix.complete <- genetic.matrix.complete[, 
                                                   order(colnames(
                                                     genetic.matrix.complete))]
morpho.matrix.complete <- morpho.matrix.complete[, 
                                                 order(colnames(
                                                   morpho.matrix.complete))]

# Prepare for analysis by preserving row names
genetic.matrix.complete.rn <- rownames_to_column(genetic.matrix.complete, 
                                                 var = "Sample")

morpho.matrix.complete.rn <- rownames_to_column(morpho.matrix.complete, 
                                                var = "Sample")

## 7.3 JACCARD SIMILARITY ANALYSIS -----

# Jaccard Index calculation function
calculate_jaccard <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union <- sum(x == 1 | y == 1)
  return(intersection / union)
}

# Calculate Jaccard Index for each sample
overlap <- sapply(1:nrow(genetic.matrix.complete.rn), function(i) {
  calculate_jaccard(genetic.matrix.complete.rn[i, -1], 
                    morpho.matrix.complete.rn[i, -1])
})

# Compile and display results
result_jaccard <- data.frame(
  Sample = genetic.matrix.complete.rn$Sample,
  Jaccard_Index = overlap
)

# Show results and summary statistics
print(result_jaccard)
mean(result_jaccard$Jaccard_Index)
sd(result_jaccard$Jaccard_Index)

# 8 DETECTION SUCCESS ----
genetic.detection <- gen.seeds.seq %>%
  select(3, 12)

genetic.detection.unique <- genetic.detection %>%
  distinct(sample, plant_genus)

genetic.detection.unique$genetic <- TRUE

morph.detection <- seeds.morph.seq %>%
  select(3,8)

morph.detection.unique <- morph.detection %>%
  distinct(sample, plant_genus)

morph.detection.unique$morph <- TRUE

# Join on sample_id and genus
detection.df <- full_join(genetic.detection.unique, 
                          morph.detection.unique, 
                          by = c("sample", "plant_genus")) %>%
                mutate(morph = ifelse(is.na(morph), FALSE, morph),
                        genetic = ifelse(is.na(genetic), FALSE, genetic)
  )

detection.summary <- detection.df %>%
  group_by(sample) %>%
  summarise(
    genetic_count = sum(genetic),
    morph_count = sum(morph),
    overlap_count = sum(genetic & morph),
    total_detected = n_distinct(plant_genus[genetic | morph]),
    method_with_more = case_when(
      genetic_count > morph_count ~ "Genetic",
      morph_count > genetic_count ~ "Morphology",
      genetic_count == morph_count ~ "Equal"
    ),
    overlap_type = case_when(
      overlap_count == 0 ~ "Absent",
      overlap_count == total_detected ~ "Total",
      overlap_count > 0 & overlap_count < total_detected ~ "Partial"
    )
  )

## 8.1 DETECTION SUMMARY COMPARISON-----

detection.summary$comparison <- ifelse(
  detection.summary$genetic_count > detection.summary$morph_count, 
  "Genetic > Morph",
  ifelse(
    detection.summary$genetic_count < detection.summary$morph_count, 
    "Morph > Genetic", 
    "Equal"
  )
)
# Create the two-way table
detection.success.overlap <- table(
  Overlap = detection.summary$overlap_type,
  Comparison = detection.summary$comparison
)

# To make it more readable, you might want to reorder the factors
detection.summary$overlap_type <- factor(
  detection.summary$overlap_type, 
  levels = c("Total", "Partial", "Absent")
)

detection.summary$comparison <- factor(
  detection.summary$comparison,
  levels = c("Equal", "Genetic > Morph", "Morph > Genetic")
)

# Now create the table with ordered factors
detection.success.overlap <- with(detection.summary, 
                                  table(overlap_type, comparison))

# Print the table
detection.success.overlap

## 8.2 LOLLIPOP GRAPH ------

order.samples.lolli <- c("297", "195",
                   "477", "027",
                   "290", 
                   "244", "087", "410", "481", "198", "199", "201", "482",
                   "255", "212", "064", "243", "242", "241", "181",
                   "384", "271", "132",
                   "253", "269",
                   "246", "103", "321", "447",
                   "113", "381",
                   "412",
                   "257", "287", "136", "110",
                   "099", "286",
                   "437", "468", "322",
                   "419", "073", "432",
                   "272", "053", "043", "070",
                   "084", "052", "050", "040"
)

detection.summary$detection_difference_count <- detection.summary$genetic_count - detection.summary$morph_count

# Reorder the dataframe
detection.summary$sample <- 
  factor(detection.summary$sample, levels = order.samples.lolli) # Apply the order

detection.difference.ordered <- detection.summary[order(
                                                  detection.summary$sample), ] 
# Reorder the rows

#visualization

ggplot(detection.difference.ordered, 
       aes(x = detection_difference_count, y = sample)) +
       geom_segment(aes(x = 0, xend = detection_difference_count, yend = sample), 
                    color = "grey") +
       geom_point(size = 3) +
       geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
       labs(
            title = 
            "Difference Between Genetic and Morphological Methods for Each Sample",
            x = "Difference (Genetic - Morphological)", y = "Sample"
  ) +
       theme_minimal() +
       theme(
            panel.grid = element_blank(),
            axis.text.y = element_text(size = 8, margin = margin(t = 10, b = 10))
  ) +
       xlim(-5, 5) +
       scale_y_discrete(expand = expansion(mult = c(0.005, 0.005)))



ggsave(
  filename = here::here("analyses", "figures&results", "lollipop_graph_COLOR.svg"), # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 8,                      # Width in inches
  height = 8,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)
## 8.3 RECALL -----

detection.summary$recall_genetic <- detection.summary$overlap_count / detection.summary$morph_count

mean(detection.summary$recall_genetic)

detection.summary$recall_morph <- detection.summary$overlap_count / detection.summary$genetic_count

mean(detection.summary$recall_morph)


## 8.4 VENN DIAGRAM -----


genetic.taxa <- df.genetic %>%
  distinct(Taxa)


seed.taxa <- df.seed %>%
  distinct(Taxa)

seed.taxa <- as.character(seed.taxa$Taxa)
genetic.taxa <- as.character(genetic.taxa$Taxa)


venn.plot <- venn.diagram(
  x = list(Morphological = seed.taxa,
           Genetic = genetic.taxa),
  category.names = c("Morphological Data", "Genetic Data"),
  filename = NULL,
  output = TRUE,
  col = "black",
  fill = c("green", "darkgreen"),
  alpha = 0.5,
  label.col = c("white", "white", "black"),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.col = c("black", "black")
)

grid.newpage()
grid.draw(venn.plot)

only.genetic <- setdiff(genetic.taxa, seed.taxa)
only.morph <- setdiff(seed.taxa, genetic.taxa)
shared <- intersect(genetic.taxa, seed.taxa)

file.remove(list.files(pattern = "*.log$"))
