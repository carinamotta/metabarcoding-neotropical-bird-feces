#reorder samples to try to make the lines easier to read

# gen_mat <- gen_mat[c("064", "181", "241", "242", "243", "253", "255",
#                      "269", "271", "468", "132", "136", "212", "246",
#                      "287", "384", "447", "482", "040", "050", "052",
#                      "084", "432", "043", "053", "070", "272",
#                      "113", "381", "321", "103", "099", "286",
#                      "322", "420", "419", "073", "110", 
#                      "199", "201", "244", "268", "290",
#                      "027", "477", "087", "410", "195", 
#                      "198", "257", "481", "297") ,]

gen_mat_2 <- gen_mat[c("040", "050", "052", "084",
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
                       "195",   "297") ,]





# gen_mat <- gen_mat[, c("Miconia", "Casearia", "Psychotria", 
#                        "Solanum", "Serjania/Paullinia",
#                        "Piper", "Trichilia", "Cestrum",
#                        "Cordia", "Alchornea", "Allophylus",
#                        "Psidium", "Morus", "Myrcia", "Nicotiana",
#                        "Cecropia", "Albizia", "Rubus", "Eugenia",
#                        "Syzygium", "Stylogyne", "Terminalia", 
#                        "Ocotea", "Ruprechtia", "Gouania", "Dendropanax",
#                        "Aspidosperma", "Maclura", "Banisteriopsis",
#                        "Senna", "Myriopus", "Tilesia", "Trema", 
#                        "Siparuna", "Lippia/Lantana", "Chamissoa", 
#                        "Centrolobium", "Amphilophium", "Citharexylum",
#                        "Guarea", "Syagrus")]

gen_mat_2 <- gen_mat_2[, c("Casearia","Ocotea", "Terminalia",
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
                           "Morus", "Rubus",  "Siparuna")]

plotweb(gen_mat_2,
        method = "cca",  # Layout method, "normal" is default
        col.high = "red",  # Color for seeds (higher level)
        col.low = "steelblue",  # Color for birds (lower level)
        text.rot = 90,  # Rotate seed labels for readability
        y.lim = c(-1, 2))  # Adjust y-axis limits if needed

#by morphospecies 
# seed_mat_morpho <- sem_morph_seq %>%
#   group_by(amostra, genus_morpho) %>%  # Group by samples and genera
#   summarise(occurrences = n(), .groups = 'drop') %>%  
#   # Count occurrences for each sample-genus pair
#   pivot_wider(names_from = amostra, values_from = occurrences,
#               values_fill = 0) %>%  # Wide format: genera as columns
#   column_to_rownames(var = "genus_morpho")


#seed_mat_morpho <- seed_mat_morpho[, c("027", "477", "195", "437",
#                                        "412", "043", "053", "052",
#                                        "385",
#                                        "040", "084", "432",
#                                        "109", "104", "050",
#                                        "070", "073", "075",
#                                        "110", "113",
#                                        "381", "384", "458",
#                                        "447", "199",
#                                        "201", "269", "064", "181",
#                                        "242", "243",
#                                        "238", "252",
#                                        "253", "255",
#                                        "271", "468", "099", "287",
#                                        "212", "321", "132", "198",
#                                        "103", "136", "246", "286",
#                                        "482", "322", "257", "420",
#                                        "481", "419", "400", "087",
#                                        "410", "244", "268", "290",
#                                        "297", "365"
#                                        )]

# seed_mat_morpho <- seed_mat_morpho[c("Morus", "Rubus_1","Pera",
#                                      "Rubus_2", "Struthanthus",
#                                      "Trichilia",
#                                      "Psychotria", 
#                                      "Casearia", "unknown_3",
#                                      "Cestrum_1", "Alchornea", 
#                                      "Cecropia", 
#                                      "Cestrum_2", "unknown_4",
#                                      "unknown_2",
#                                      "Allophylus",
#                                      "Miconia_1", "Myrtaceae_3",
#                                      "Myrtaceae_2",
#                                      "Solanum_2",
#                                      "Piper_2",
#                                      "Miconia_3", "Nicotiana",
#                                      "Myriopus", 'Maclura', 
#                                      "Lantanta", "Solanum_1",
#                                      "Miconia_2",
#                                      "Tilesia", "Stylogyne",
#                                      "Piper_1", "Amaranthaceae",
#                                      "Myrcia",
#                                      "Fabaceae_1",
#                                      "Miconia_4", "Psidium_1",
#                                      "unknown_1","Psidium_2",
#                                      "Siparuna", "unknown_5"), ]

# plotweb(seed_mat_morpho,
#         method = "normal",  # Layout method, "normal" is default
#         col.high = "steelblue",  # Color for seeds (higher level)
#         col.low = "darkgreen",  # Color for birds (lower level)
#         text.rot = 90,  # Rotate seed labels for readability
#         y.lim = c(-1, 2))  # Adjust y-axis limits if needed



#reorder samples to make lines easier to read
seed_mat_genus <- seed_mat_genus[, c("040", "050", "052", "084",
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
)]



seed_mat_genus <- seed_mat_genus[c( "Casearia","Trichilia","unknown_3",
                                    
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
                                    "Siparuna"), ]



plotweb(seed_mat_genus,
        method = "cca",  # Layout method, "normal" is default
        col.high = "steelblue",  # Color for seeds (higher level)
        col.low = "darkblue",  # Color for birds (lower level)
        text.rot = 90,  # Rotate seed labels for readability
        y.lim = c(-1, 2))  # Adjust y-axis limits if needed

#by genus, switched columns and rows to match with genetic matrix 
seed_mat_genus <- sem_morph_seq %>%
  group_by(amostra, genus) %>%  # Group by samples and genera
  summarise(occurrences = n(), .groups = 'drop') %>%  # Count occurrences for each sample-genus pair
  pivot_wider(names_from = amostra, values_from = occurrences,
              values_fill = 0) %>%  # Wide format: genera as columns
  column_to_rownames(var = "genus")
