
bipartite_data <- data.frame(
  Genera_Left = c("GenusA", "GenusA", "GenusB", "GenusC", "GenusD"),
  Samples = c("Sample1", "Sample1", "Sample2", "Sample2", "Sample3"),
  Genera_Right = c("GenusX", NA, "GenusZ", "GenusX", "GenusY")  # NA for unmatched
)

install.packages("ggalluvial")
library(ggalluvial)

ggplot(bipartite_data, aes(axis1 = Genera_Left, axis2 = Samples, axis3 = Genera_Right)) +
  geom_alluvium(aes(fill = Samples), width = 0.2, knot.pos = 0.5) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera_Left", "Samples", "Genera_Right"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  labs(
    title = "Genera-Sample-Genera Alluvial Diagram (Presence/Absence)",
    x = "Groups",
    y = "Presence"
  )



# Example presence/absence data
bipartite_data <- data.frame(
  x = c("GenusA", "GenusA", "GenusB", "Sample1", "Sample1", "Sample2", "GenusX", "GenusZ"),
  next_x = c("Sample1", "Sample1", "Sample2", "GenusX", "GenusY", "GenusZ", "GenusX", "GenusZ"),
  node = c(1, 1, 1, 1, 1, 1, 1, 1)  # Presence/absence represented as 1
)




df_seed <- subset(seed_complete, Weight >= 1)


desired_order_seed <- c( "Miconia", "Casearia", "Piper", 
                         "Solanum", "Psychotria", "Myriopus",
                         "Trichilia", "Cestrum", "Cecropia",
                         "Rubus", "Allophylus", "Fabaceae_1",
                         "Psidium", "unknown_4", "Stylogyne",
                         "Morus", "Myrcia", "Nicotiana", "Alchornea",
                         "Pera", "Struthanthus", "Lantanta", "Maclura", 
                         "Myrtaceae_2", "unknown_1", "Tilesia", "Siparuna",
                         "unknown_2", "Amaranthaceae", "unknown_3")



df_seed$Genus <- factor(df_seed$Genus, levels = desired_order_seed)


# Create a custom palette with 20 tones of purple
purple_palette <- colorRampPalette(c("#8A2BE2", "Purple"))(40)

# View the palette
print(purple_palette)

# Plot to visualize the colors
barplot(rep(1, 40), col = purple_palette, border = NA, space = 0, main = "20 Tones of Purple")

# Create a custom palette with 20 tones of violet
violet_palette <- colorRampPalette(c("#8A2BE2", "#EE82EE"))(20)  # From BlueViolet to Violet

# View the palette
print(violet_palette)

# Plot to visualize the colors
barplot(rep(1, 20), col = violet_palette, border = NA, space = 0, main = "20 Tones of Violet")


# Define a gradient palette for hot pink to red
hot_pink_red_palette <- colorRampPalette(c("#FF69B4", "#FF0000"))  # Hot pink to red

# Generate a palette with a specific number of colors (e.g., 20)
gradient_colors <- hot_pink_red_palette(40)

print(gradient_colors)


# Define a gradient palette for bright blue colors
blue_gradient_palette <- colorRampPalette(c("#1E90FF", "#005A9C"))  # Sky Blue to Royal Blue

# Generate a palette with a specific number of colors (e.g., 20)
gradient_blue_colors <- blue_gradient_palette(20)

# Visualize the palette
barplot(rep(1, 20), col = gradient_blue_colors, border = NA, space = 0, main = "Bright Blue Gradient")

print(gradient_blue_colors)

# Visualize the palette
barplot(rep(1, 40), col = gradient_colors, border = NA, space = 0, main = "Hot Pink to Red Gradient")



# Install and load the RColorBrewer package if you haven't already
# install.packages("RColorBrewer")
library(RColorBrewer)

# Use a color-blind friendly palette with 20 unique colors
color_blind_friendly_palette <- brewer.pal(12, "Set3")  # This palette has 12 colors
# If you want 20 colors, you could combine multiple color sets
# For example, using 2 sets of "Set3" and adding "Dark2" for variety
extended_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))

# Create a vector of your original 8 colors
colors <- c("#192609", "#75774E", "#F19E14", "#6A0207", "#71B481", "#3B5708", "#5B7C91", "#9F4147")

# Use colorRampPalette to generate 20 interpolated colors
color_palette <- colorRampPalette(colors)(20)

# Print the new palette
color_palette

barplot(rep(1, 20), col = color_palette, border = NA, space = 0)

print(color_palette)


list(extended_palette
     )
# Visualize the palette
barplot(rep(1, length(extended_palette)), col = extended_palette, border = NA, space = 0, main = "Color-Blind Friendly Palette")

ggplot(df_seed,
       aes(axis1 = Genus, axis2 = Sample, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.2) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Genera-Sample Alluvial Diagram",
    x = "Groups",
    y = "Presence"
  )

ggsave(
  filename = "morph_graph.svg", # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 4,                      # Width in inches
  height = 16,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)

####-------

df_seed <- subset(seed_complete, Weight >= 1)

desired_order_seed_samples_2 <-  c("040", "050", "052", "084",
                                     "070","043", "053", "272", 
                                     "432","420", "073","419", 
                                     "322","468", "437",
                                     
                                     "286","099",
                                     "110", "136", "287", "257",
                                     "412",
                                     "381","113",
                                     "447", "321", "103","246",
                                     "269", "253", 
                                     "132", "271",  "384", 
                                     "181","241", "242", "243", "064", "212", "255",
                                     "482",
                                     "201", "199","198","481", "410","087", "244",
                                     "268", "290", 
                                     "027", "477", 
                                     "195",   "297"
)

desired_order_seed <- c( "Casearia","Trichilia","unknown_3",
                           
                           "Piper", "unknown_4", "Cestrum",
                           "Chamissoa",
                           "Amaranthaceae",
                           "Nicotiana",
                           "Tilesia",
                           "Alchornea","Cecropia", 
                           "unknown_2",
                           "Psychotria",
                           "Myriopus",
                           "Lantanta",
                           "Solanum", "Miconia",
                           
                           
                           "Allophylus", "Maclura", "Myrtaceae_2",
                           
                           
                           "Stylogyne", 
                           "Fabaceae_1",
                           "Myrcia", "unknown_1","Psidium","Morus",
                           
                           
                           
                           
                           "Rubus","Struthanthus","Pera",
                           "Siparuna"
)


# Reorder the Genera column
df_seed$Sample <- factor(df_seed$Sample, levels = desired_order_seed_samples_2)

df_seed$Genus <- factor(df_seed$Genus, levels = desired_order_seed_2)



ggplot(df_seed,
       aes(axis1 = Genus, axis2 = Sample, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.1, alpha = 0.7) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Morus" = "#192609",  # Dark Green
      "Casearia" = "#577882",  # Slate Blue
      "Trichilia" = "#5C613B",  # Muted Green
      "Miconia" = "#9F4147",   # Deep Red
      "Cestrum" = "#AF8932",  # Amber
      "Myrcia" = "#DD971D",  # Golden Yellow
      "Nicotiana" = "#D47D11",  # Dark Orange
      "Solanum" = "#A2430C",  # Rust Red
      "Alchornea" = "#710A07",  # Burgundy
      "Cecropia" = "#6C3A2D",  # Warm Brown
      "Psychotria" = "#6E7B5A",  # Sage Green
      "Rubus" = "#6EAF7A",  # Soft Green
      "Allophylus" = "#5A8C4E",  # Grass Green
      "Maclura" = "#466A21",  # Forest Green
      "Psidium" = "#405C1D",  # Olive Green
      "Piper" = "#4B6A50",  # Moss Green
      "Myriopus" = "#3A4322",  # Olive Green
      "Tilesia" = "#6C6C7D",  # Grayish Blue
      "Siparuna" = "#855662",  # Muted Purple
      "Stylogyne" = "#827B47",  # Olive Brown
      "Fabaceae_1" = "black",  # Bright Blue
      "unknown_4" = "black",  # Bright Blue
      "Pera" = "black",  # Bright Blue
      "Struthanthus" = "black",  # Bright Blue
      "Lantanta" = "black",  # Bright Blue
      "Myrtaceae_2" = "black",  # Bright Blue
      "unknown_1" = "black",  # Bright Blue
      "unknown_2" = "black",  # Bright Blue
      "Amaranthaceae" = "black",  # Bright Blue
      "unknown_3" = "black"   # Bright Blue
    )
  ) 

# 
# ggplot(df_seed,
#        aes(axis1 = Genus, axis2 = Sample, y = Weight)) +
#   geom_alluvium(aes(fill = Genus), width = 0.2) +
#   geom_stratum(width = 0.2, fill = "gray", color = "black") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_manual(
#     values = c(
#       "Morus" = "#1B9E77",  # Green
#       "Casearia" = "#61864B",  # Olive Green
#       "Trichilia" = "#A66F20",  # Warm Yellow
#       "Miconia" = "#CE6014",  # Orange
#       "Cestrum" = "#A96755",  # Brownish Red
#       "Myrcia" = "#846D97",  # Lavender
#       "Nicotiana" = "#8D61AA",  # Purple-Lavender
#       "Solanum" = "#B7469B",  # Purple-Pink
#       "Alchornea" = "#E12C8C",  # Bright Pink
#       "Cecropia" = "#BE5067",  # Deep Pink
#       "Psychotria" = "#8E7E40",  # Olive Brown
#       "Rubus" = "#6CA61C",  # Olive Green
#       "Allophylus" = "#9BA812",  # Olive Yellow
#       "Maclura" = "#CBA907",  # Gold
#       "Psidium" = "#DBA206",  # Dark Yellow
#       "Piper" = "#C48F10",  # Warm Yellow-Orange
#       "Myriopus" = "#AC7B1A",  # Yellow-Brown
#       "Tilesia" = "#957130",  # Mustard Yellow
#       "Siparuna" = "#7D6B4B",  # Brownish Green
#       "Stylogyne" = "#666666",  # Gray
#       
#       "Fabaceae_1" = "darkgray",  # Bright Blue
#       "unknown_4" = "#147EDF",  # Bright Blue
#       "Pera" = "#127CDA",  # Bright Blue
#       "Struthanthus" = "#1179D5",  # Bright Blue
#       "Lantanta" = "#0F76D0",  # Bright Blue
#       "Myrtaceae_2" = "#0E73CA",  # Bright Blue
#       "unknown_1" = "#0C70C5",  # Bright Blue
#       "unknown_2" = "#0B6DC0",  # Bright Blue
#       "Amaranthaceae" = "#096BBB",  # Bright Blue
#       "unknown_3" = "#0768B6"   # Bright Blue
#     )
#   ) +
#   labs(
#     title = "Genera-Sample Alluvial Diagram",
#     x = "Groups",
#     y = "Presence"
#   )


ggsave(
  filename = "morph_graph_3.svg", # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 4,                      # Width in inches
  height = 16,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)

###---------

df_gen <- subset(gen_complete, Weight >= 1)



ggplot(df_gen,
       aes(axis1 = Sample, axis2 = Genus, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.2) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Genera-Sample Alluvial Diagram",
    x = "Groups",
    y = "Presence"
  )

ggsave(
  filename = "gen_graph.svg", # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 4,                      # Width in inches
  height = 16,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)



desired_order_gen <- c( "Miconia", "Casearia", "Psychotria", 
                        "Solanum", "Serjania/Paullinia",
                                                "Piper", "Trichilia", "Cestrum",
                                                "Cordia", "Alchornea", "Allophylus",
                                                "Psidium", "Morus", "Myrcia", "Nicotiana",
                                                "Cecropia", "Albizia", "Rubus", "Eugenia",
                                                "Syzygium", "Stylogyne", "Terminalia", 
                                                "Ocotea", "Ruprechtia", "Gouania", "Dendropanax",
                                                "Aspidosperma", "Maclura", "Banisteriopsis",
                                                "Senna", "Myriopus", "Tilesia", "Trema", 
                                                "Siparuna", "Lippia/Lantana", "Chamissoa", 
                                                "Centrolobium", "Amphilophium", "Citharexylum",
                                                "Guarea", "Syagrus"
)


# Reorder the Genera column
#df_gen$Sample <- factor(df_gen$Sample, levels = desired_order)

df_gen$Genus <- factor(df_gen$Genus, levels = desired_order_gen)


####------

df_gen <- subset(gen_complete, Weight >= 1)

desired_order_samples <- c("040", "050", "052", "084",
                           "070","043", "053", "272", 
                           "432","420", "073","419", 
                           "322","468", "437",
                           
                           "286","099",
                           "110", "136", "287", "257",
                           "412",
                           "381","113",
                           "447", "321", "103","246",
                           "269", "253", 
                           "132", "271",  "384", 
                           "181","241", "242", "243", "064", "212", "255",
                           "482",
                           "201", "199","198","481", "410","087", "244",
                           "268", "290", 
                           "027", "477", 
                           "195",   "297"
)


desired_order_gen <- c( "Casearia","Ocotea", "Terminalia",
                        "Trichilia", "Cestrum",
                        "Chamissoa","Piper",
                        
                        "Cecropia", 
                        "Serjania/Paullinia","Centrolobium",
                        "Myriopus", "Tilesia",
                        "Ruprechtia", 
                        "Gouania","Citharexylum",
                        "Solanum","Syzygium","Banisteriopsis",
                        "Albizia", "Aspidosperma",
                        "Lippia/Lantana","Psychotria",
                        
                        "Alchornea",
                        "Nicotiana", "Amphilophium",
                        
                        "Cordia",
                        "Miconia","Trema", "Senna",
                        
                        
                        
                        
                        
                        "Maclura",
                        "Dendropanax",  
                        "Stylogyne","Allophylus",
                        "Eugenia",
                        "Guarea", "Syagrus", 
                        "Myrcia","Psidium",
                        "Morus", "Rubus",  "Siparuna"
)


# Reorder the Genera column
df_gen$Sample <- factor(df_gen$Sample, levels = desired_order_samples)

df_gen$Genus <- factor(df_gen$Genus, levels = desired_order_gen)

ggplot(df_gen,
       aes(axis1 = Sample, axis2 = Genus, y = Weight)) +
  geom_alluvium(aes(fill = Genus), width = 0.1, alpha = 0.7) +
  geom_stratum(width = 0.2, fill = "gray", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Genera", "Sample"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Morus" = "#192609",  # Dark Green
      "Casearia" = "#577882",  # Slate Blue
      "Trichilia" = "#5C613B",  # Muted Green
      "Miconia" = "#9F4147",   # Deep Red
      "Cestrum" = "#AF8932",  # Amber
      "Myrcia" = "#DD971D",  # Golden Yellow
      "Nicotiana" = "#D47D11",  # Dark Orange
      "Solanum" = "#A2430C",  # Rust Red
      "Alchornea" = "#710A07",  # Burgundy
      "Cecropia" = "#6C3A2D",  # Warm Brown
      "Psychotria" = "#6E7B5A",  # Sage Green
      "Rubus" = "#6EAF7A",  # Soft Green
      "Allophylus" = "#5A8C4E",  # Grass Green
      "Maclura" = "#466A21",  # Forest Green
      "Psidium" = "#405C1D",  # Olive Green
      "Piper" = "#4B6A50",  # Moss Green
      "Myriopus" = "#3A4322",  # Olive Green
      "Tilesia" = "#6C6C7D",  # Grayish Blue
      "Siparuna" = "#855662",  # Muted Purple
      "Stylogyne" = "#827B47",  # Olive Brown
      "Ocotea" = "#000000",  # Black
      "Terminalia" = "#000000",  # Black
      "Chamissoa" = "#000000",  # Black
      "Serjania/Paullinia" = "#000000",  # Black
      "Centrolobium" = "#000000",  # Black
      "Ruprechtia" = "#000000",  # Black
      "Gouania" = "#000000",  # Black
      "Citharexylum" = "#000000",  # Black
      "Syzygium" = "#000000",  # Black
      "Banisteriopsis" = "#000000",  # Black
      "Albizia" = "#000000",  # Black
      "Aspidosperma" = "#000000",  # Black
      "Lippia/Lantana" = "#000000",  # Black
      "Amphilophium" = "#000000",  # Black
      "Cordia" = "#000000",  # Black
      "Trema" = "#000000",  # Black
      "Senna" = "#000000",  # Black
      "Dendropanax" = "#000000",  # Black
      "Eugenia" = "#000000",  # Black
      "Guarea" = "#000000",  # Black
      "Syagrus" = "#000000"   # Black
      
    ))
            +
  labs(
    title = "Genera-Sample Alluvial Diagram",
    x = "Groups",
    y = "Presence"
  )

ggsave(
  filename = "gen_graph_3.svg", # File name
  plot = last_plot(),             # Optional: explicitly specify the plot
  width = 4,                      # Width in inches
  height = 16,                     # Height in inches
  dpi = 720                       # Resolution in dots per inch
)


######-----



gen_morph <- merge(x=df_gen, y=df_seed, by="Genus", all.x =T)


gen_morph <- gen_morph %>%
  filter(!is.na(Sample.y))

gen_morf_list <- gen_morph %>% distinct(Genus, .keep_all = F)

list(gen_morf_list)
