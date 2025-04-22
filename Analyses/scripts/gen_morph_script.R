## Created: 02-October-2024
## Updated: 18-April-2025

## Author: Carina Isabella Motta

## Evaluating a  DNA metabarcoding approach of bird feces to quantify and predict seed dispersal

# Data cleaning and analysis to compare two data sets: Morphological identifications of seeds and identification of plant DNA in same bird feces 

# 1 LOAD PACKAGES------------------------------------------------------------

# a vector listing package names needed 

package.list <- c("here", #makes this project the working directory
                  "tidyverse", #data cleaning
                  "dplyr", #data cleaning
                  "reshape2", 
                  "igraph",
                  "VennDiagram", #VennDiagram figure
                  "ggalluvial" 
)


#installing the packages if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()
                               [,"Package"])]

if(length(new.packages)) install.packages(new.packages)

#and loading the packages into R with a for loop
for(i in package.list){library(i, character.only = T)}

# 2 LOAD DATA---------------------------------------------------------------

## 2.1 LOAD SEED DATA----

sem <- readr::read_csv(here::here("Analyses","data",                                "seeds_28MAR2025_2.csv"))

#transform into a tibble to facilitate using packages like dyplr 
sem <- as_tibble(sem)

## 2.2 LOAD SEED MORPHOSPECIES DATA----

morpho <- readr::read_csv(here::here("Analyses", "data",                               "morpho_4APRIL2025.csv"))

#transform into a tibble to faciliate using packages like dyplr 
morpho <- as_tibble(morpho)

## 2.3 LOAD METABARCODING DATA----

gen <- readr::read_csv(here::here("Analyses","data",                                          "metabarcoding_results.csv"))

#transform into a tibble to faciliate using packages like dyplr 
gen <- as_tibble(gen)

# 3 DATA CLEANING AND PREP-------------------------------------------------------

## 3.1 SEED & SEED MORPHOLOGY DATA----

### 3.1.1 SEED DATA-----

# subset columns that we are interested in 
sem <- sem %>% 
  select(3, 4, 6, 7, 8)

### 3.1.2 MORPHOLOGY DATA-----

# subset columns that we are interested in 
morpho <- morpho %>% 
  select(1:5)

### 3.1.3 MERGE SEED AND SEED MORPHOLOGY DATA----

sem_morph <- merge(x=sem, y=morpho, by="morpho_spp", all.x =T)

##3.2 METABARCODING DATA----

# cut excess rows and columns 

#cut excess rows
gen <- gen[-(561:1246),]

#cut excess columns 
gen <- gen[,-(16:32)]

#filter to only include pools that include samples with seeds, and corresponding blanks

gen_sem <- gen %>%
  filter(pool %in% c("P_00", "P_01", "P_02", "P_03","P_04",                         "P_05", "P_06", "P_07", "P_08", "P_09",                        "P_10", "P_11", "P_12"))%>%
  filter(seeds == 1 | grepl("blank", sample))


# 4 SEED & SEED MORPHOLOGY SUMMARY STATISTICS----
## 4.1 ALL SAMPLES----


n_distinct(sem_morph$unique)

n_distinct(sem_morph$plant_family)

n_distinct(sem_morph$plant_genus)

sum(sem_morph$no_seeds)

bird.spp.morph <- sem_morph %>%
                  group_by(bird_species) %>%
                  summarize(n = n_distinct(sample))

## 4.2 SUCESSFULLY SEQUENCED SAMPLES----
gen_sem_seq <- gen_sem %>% 
               filter(final_fate == 1)

#cross reference seed morphology data set with sequencing to only include 
#the samples that were successfully sequenced. That way we can create our 
#interaction webs using only data that is present in both data sets

sem_morph_seq <- sem_morph %>%
  filter(sample %in% gen_sem_seq$sample)

#check to see if was correctly filtered, should have 52 samples now 

n_distinct(sem_morph_seq$sample) # 52, correctly filtered 

#check the number of morphospecies present after filtering out the samples 
#that didn't successfully sequence, should be 36

n_distinct(sem_morph_seq$unique)

#create a table of how many samples of each species was successfully 
#sequenced 

spp_gen_suc <- sem_morph_seq %>%
  group_by(bird_species) %>%
  summarize(n = n_distinct(sample))

n_distinct(sem_morph_seq$unique)

n_distinct(sem_morph_seq$plant_family)

n_distinct(sem_morph_seq$plant_genus)

sum(sem_morph_seq$no_seeds)


### 4.2.1----

# 5 METABARCODING SUMMARY STATISTICS----
## 5.1 FATE SUMMARIES----

# final fates of each OTU
# 0 means it didn't make the read threshold,
# 1 means is was successful and kept
# 2 means it didn't make the match threshold
# 3 means it wasn't successfully sequenced 
# 4 means it was excluded due to suspected contamination 
# 5 means that it was a blank sample that generated an OTU

sample_fate <- gen_sem %>%
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
#       while individual OTUs were deleted due to insufficient match or 
#       contamination, those two were not motives for excluding samples
sample_fate_summary <- sample_fate %>%
  group_by(final_fate) %>%
  summarise(total_samples = n())


OTU_fate <- gen_sem %>%
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

OTU_fate_summary <- OTU_fate %>%
  group_by(final_fate) %>%
  summarise(total_OTUs = n_distinct(OTU), n_reads = sum(no_reads_OTU))

## 5.2 READS & OTUs before quality filtering ---- 
t_reads_OTU <- gen_sem %>%
  filter(!grepl("blank", sample),!is.na(OTU))  %>%
  group_by(OTU) %>%
  summarise(t_reads = sum(no_reads_OTU))

#the sum of reads
sum(t_reads_OTU$t_reads)

#the mean number of reads per OTU
mean(t_reads_OTU$t_reads)

#the standard deviation of reads
sd(t_reads_OTU$t_reads)

## 5.3 READS & OTUs after quality filtering ----- 
t_reads_OTU_filtered <- gen_sem %>%
  filter(final_fate == 1)  %>%
  group_by(OTU) %>%
  summarise(t_reads = sum(no_reads_OTU))

#the sum of reads
sum(t_reads_OTU_filtered$t_reads)

#the mean number of reads per OTU
mean(t_reads_OTU_filtered$t_reads)

#the standard deviation of reads
sd(t_reads_OTU_filtered$t_reads)

#per dropping sample
t_reads_sample_filtered <- gen_sem %>%
  filter(final_fate == 1)  %>%
  group_by(sample) %>%
  summarise(t_reads_sample = sum(no_reads_OTU))

## 5.3 PLANT GENERA----

n_plant_fam <- gen_sem %>%
  filter(!grepl("blank", sample), !is.na(OTU)) %>%
  filter(final_fate == 1)  %>%
  summarise(n_plant_fam = n_distinct(plant_family), 
            n_plant_genera = n_distinct(plant_genus))

n_plant_fam$n_plant_fam

# 6 TRIPARTITE VISUALIZATION ----
## 6.1 CREATE INTERACTION MATRIX WITH GENETIC DATA-----

# Create genetic matrix
gen_mat <- gen_sem_seq %>%
  group_by(sample, plant_genus) %>%
  summarise(occurrences = n(), .groups = 'drop') %>%
  pivot_wider(names_from = plant_genus, 
              values_from = occurrences, 
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

# Convert to presence/absence (1/0)
gen_mat[gen_mat > 1] <- 1

# Prepare for visualization


# Add rownames as column for melting
gen_mat_with_rowname <- rownames_to_column(gen_mat, var = "Sample")

# Transform to long format
gen_list <- melt(gen_mat_with_rowname)
colnames(gen_list) <- c("Sample", "Genus", "Weight")

# Ensure presence/absence (1/0)
gen_list$Weight[gen_list$Weight > 1] <- 1

# Filter to only positive detections
df_gen <- subset(gen_list, Weight >= 1)

# Define factor orders

# Sample order
desired_order_samples <- c(
  "040", "050", "052", "084", "070", "043", "053", "272", 
  "432", "073", "419", "322", "468", "437", "286", "099",
  "110", "136", "287", "257", "412", "381", "113", "447", 
  "321", "103", "246", "269", "253", "132", "271", "384", 
  "181", "241", "242", "243", "064", "212", "255", "482",
  "201", "199", "198", "481", "410", "087", "244", "290", 
  "027", "477", "195", "297"
)

# Genus order
desired_order_gen <- c(
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
df_gen$Sample <- factor(df_gen$Sample, levels = desired_order_samples)
df_gen$Genus <- factor(df_gen$Genus, levels = desired_order_gen)

# Create alluvial plot

gen_plot <- ggplot(df_gen,
                   aes(axis1 = Sample, axis2 = Genus, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.1, alpha = 0.7) +
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

gen_plot

# Save outputs

ggsave(
  filename = "gen_graph.svg",
  plot = gen_plot,
  width = 4,
  height = 16,
  dpi = 720
)

## 6.2 CREATE INTERACTION MATRIX WITH MORPHOLOGY DATA-----

# Create morphology matrix

seed_mat <- sem_morph_seq %>%
  group_by(sample, plant_genus) %>%
  summarise(occurrences = n(), .groups = 'drop') %>%
  pivot_wider(names_from = plant_genus, 
              values_from = occurrences,
              values_fill = 0) %>%
  column_to_rownames(var = "sample")

# Prepare for visualization

# Add rownames as column for melting
seed_mat_with_rowname <- rownames_to_column(seed_mat, var = "Sample")

# Transform to long format
seed_list <- melt(seed_mat_with_rowname)
colnames(seed_list) <- c("Sample", "Genus", "Weight")

# Ensure presence/absence (1/0)
seed_list$Weight[seed_list$Weight > 1] <- 1

# Save intermediate matrix
#write.csv(seed_mat, here::here("figures&results", "seed_mat.csv"))

# Filter to only positive detections
df_seed <- subset(seed_list, Weight >= 1)

# Define factor orders

# Genus order for morphology data
desired_order_seed <- c(
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
df_seed$Sample <- factor(df_seed$Sample, levels = desired_order_samples)
df_seed$Genus <- factor(df_seed$Genus, levels = desired_order_seed)


# Create alluvial plot

morph_plot <- ggplot(df_seed,
                     aes(axis1 = Genus, axis2 = Sample, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.1, alpha = 0.7) +
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

morph_plot

# save outputs 
ggsave(
  filename = "morph_graph.svg",
  plot = morph_plot,
  width = 4,
  height = 16,
  dpi = 720
)

write.csv(df_seed, here::here("figures&results", "df_seed.csv"))
# 7 JACCARD INDEX ----
## 7.1 CREATE MATRICES OF THE TWO DATASETS WITH MATCHING DIMENSIONS----

# Convert genetic matrix to long format
gen_mat_with_rowname <- rownames_to_column(gen_mat, var = "Sample")
gen_list <- melt(gen_mat_with_rowname)
colnames(gen_list) <- c("Sample", "Genus", "Weight")

# Convert to presence/absence (1/0)
gen_list$Weight[gen_list$Weight > 1] <- 1

# Convert morphology matrix to long format
seed_mat_with_rowname <- rownames_to_column(seed_mat, var = "Sample")
seed_list <- melt(seed_mat_with_rowname)
colnames(seed_list) <- c("Sample", "Genus", "Weight")

# Convert to presence/absence (1/0)
seed_list$Weight[seed_list$Weight > 1] <- 1

# Show ordered summary of genetic detections
gen_list_summary <- gen_list %>%
  group_by(Genus) %>%
  summarise(n = sum(Weight))

ordered <- gen_list_summary[order(-gen_list_summary$n), ]
print(ordered, n = 41)

## 7.2 ALIGN MATRIX DIMENSIONS BY HANDLING MISSING DATA ----

# FIND AND HANDLE DATA MISSING FROM MORPHOLOGY DATASET
missing_from_morph <- gen_list %>% 
  anti_join(seed_list, by = "Genus")
missing_from_morph$Weight <- 0
seed_complete <- bind_rows(seed_list, missing_from_morph)

# FIND AND HANDLE DATA MISSING FROM GENETIC DATASET
missing_from_gen <- seed_list %>% 
  anti_join(gen_list, by = "Genus")
missing_from_gen$Weight <- 0
gen_complete <- bind_rows(gen_list, missing_from_gen)

# RECONSTRUCT MATRICES WITH COMPLETE DIMENSIONS
genetic_matrix <- gen_complete %>%
  pivot_wider(names_from = Genus, values_from = Weight) %>%
  column_to_rownames(var = "Sample")

morpho_matrix <- seed_complete %>%
  pivot_wider(names_from = Genus, values_from = Weight) %>%
  column_to_rownames(var = "Sample")

# Standardize column order
genetic_matrix <- genetic_matrix[, order(colnames(genetic_matrix))]
morpho_matrix <- morpho_matrix[, order(colnames(morpho_matrix))]

# Prepare for analysis by preserving row names
genetic_matrix_w_rowname <- rownames_to_column(genetic_matrix, var = "Sample")
morpho_matrix_w_rowname <- rownames_to_column(morpho_matrix, var = "Sample")

## 7.3 JACCARD SIMILARITY ANALYSIS -----

# Jaccard Index calculation function
calculate_jaccard <- function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union <- sum(x == 1 | y == 1)
  return(intersection / union)
}

# Calculate Jaccard Index for each sample
overlap <- sapply(1:nrow(genetic_matrix_w_rowname), function(i) {
  calculate_jaccard(genetic_matrix_w_rowname[i, -1], 
                    morpho_matrix_w_rowname[i, -1])
})

# Compile and display results
result_jaccard <- data.frame(
  Sample = genetic_matrix_w_rowname$Sample,
  Jaccard_Index = overlap
)

# Show results and summary statistics
print(result_jaccard)
mean(result_jaccard$Jaccard_Index)
sd(result_jaccard$Jaccard_Index)

# 8 DETECTION SUCCESS ----
genetic_data <- gen_sem_seq %>%
  select(3, 12)

unique_detections_genetic <- genetic_data %>%
  distinct(sample, plant_genus)

unique_detections_genetic$genetic <- TRUE

morph_data <- sem_morph_seq %>%
  select(3,8)

unique_detections_morph <- morph_data %>%
  distinct(sample, plant_genus)

unique_detections_morph$morph <- TRUE

# Join on sample_id and genus
detection_df <- full_join(unique_detections_genetic, unique_detections_morph, by = c("sample", "plant_genus")) %>%
  mutate(
    morph = ifelse(is.na(morph), FALSE, morph),
    genetic = ifelse(is.na(genetic), FALSE, genetic)
  )

detection_summary <- detection_df %>%
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

detection_summary$detection_success_morph <- detection_summary$morph_count / detection_summary$total_detected

mean(detection_summary$detection_success_morph)
sd(detection_summary$detection_success_morph)

detection_summary$detection_success_genetic <- detection_summary$genetic_count / detection_summary$total_detected

mean(detection_summary$detection_success_genetic)
sd(detection_summary$detection_success_genetic)

mean(detection_summary$overlap_type == "Total")



detection_summary$comparison <- ifelse(
  detection_summary$genetic_count > detection_summary$morph_count, "Genetic > Morph",
  ifelse(
    detection_summary$genetic_count < detection_summary$morph_count, "Morph > Genetic", 
    "Equal"
  )
)
# Create the two-way table
result_table <- table(
  Overlap = detection_summary$overlap_type,
  Comparison = detection_summary$comparison
)

# To make it more readable, you might want to reorder the factors
detection_summary$overlap_type <- factor(
  detection_summary$overlap_type, 
  levels = c("Total", "Partial", "Absent")
)

detection_summary$comparison <- factor(
  detection_summary$comparison,
  levels = c("Equal", "Genetic > Morph", "Morph > Genetic")
)

# Now create the table with ordered factors
result_table <- with(detection_summary, table(overlap_type, comparison))

# Print the table
result_table

## 8.1 Lollipop graph visualization 

order_samples <- c("297", "195",
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

detection_summary$detection_difference <- detection_summary$detection_success_genetic - detection_summary$detection_success_morph
# Reorder the dataframe
detection_summary$sample <- 
  factor(detection_summary$sample, levels = order_samples) # Apply the order

recall_difference_df_ordered <- detection_summary[order(detection_summary$sample), ] # Reorder the rows


#visualization

ggplot(recall_difference_df_ordered, aes(x = detection_difference, y = sample)) +
  geom_segment(aes(x = 0, xend = detection_difference, y = sample, yend = sample), color = "grey") +  # Horizontal line
  geom_point(aes(x = detection_difference, y = sample), 
             size = 3) +  # Dot at the end
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Vertical line at 0
  labs(title = "Difference Between Genetic and Morphological Methods for Each Sample",
       x = "Difference (Genetic - Morphological)",
       y = "Sample") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove all gridlines
        axis.text.y = element_text(size = 8)) +
  xlim(-3, 3) +
  scale_y_discrete(expand = expansion(mult = c(0.005, 0.005))) +
  theme(axis.text.y = element_text(margin = margin(t = 10, b = 10)))


ggsave(
  filename = "lollipop_graph_COLOR.svg", # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 8,                      # Width in inches
  height = 8,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)

# 9 VENN DIAGRAM -----

n_distinct(detection_df$plant_genus)

list_gen <- gen_list %>%
  filter(Weight == 1) %>%
  distinct(Genus)

print(list_gen)


list_seeds <- seed_list %>%
  filter(Weight == 1) %>%
  distinct(Genus)

n_distinct(list_seed_onlygen$Genus)

list_seed_onlygen <- list_seeds[-c(18, 24:28), ]

venn.plot <- venn.diagram(
  x = list(Genetic = gen_list$Genus, 
           Morphological = list_seed_onlygen$Genus),
  category.names = c("Genetic Data", "Morphological Data"),
  filename = NULL,
  output = TRUE,
  col = "black",
  fill = c("darkgreen", "orange"),
  alpha = 0.5,
  label.col = c("white", "white", "black"),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontfamily = "sans",
  cat.col = c("black", "black")
)

grid.draw(venn.plot)

only_genetic <- setdiff(list_gen, list_seeds)
only_morph <- setdiff(list_seeds, list_gen)
shared <- intersect(list_gen, list_seeds)

writeClipboard(as.character(only_genetic$Genus))

writeClipboard(as.character(only_morph$Genus))

writeClipboard(as.character(shared$Genus))

