## Created: 02-October-2024
## Updated: 26-November-2024

## Carina Isabella Motta

## Evaluating a  DNA metabarcoding approach of bird feces to quantify and 
## predict seed dispersal

# Data cleaning and analysis to compare two data sets: Morphological 
# identifications of seeds and identification of plant DNA in same bird feces 

# Proposed stats and analyses (see SECTION 4)------

#   1. Number of samples (with seeds) analyzed 
#   2. Number of samples successfully sequenced 
#   3. OTU stats:
#       3a. Combined read count, total number of paired reads
#       3b. Mean paired reads / OTU
#       3c. Number of OTUs produced 
#       3d. Number of OTUs maintained 
#           3d1. % discarded due to insufficient reads 
#           3d2. % discarded due to insufficient match 
#       3e. Number of paired reads maintained / OTU
#       3f. Mean number of OTU / ind
#       3g. Mean number of OTU / spp.
#       3h. Number of plant families identified 
#       3i. Number of plant genera identified 
#       3j. Number of detections / OTU --> accumulation curve
#   4. Tripartite network analysis 
#   5. Recall or directional analysis 


#1 LOAD PACKAGES----------------------------------------------------------------
  
  # a vector listing package names needed 
  
  package.list <- c("here", #makes this project the working directory
                    "tidyverse", #data cleaning
                    "dplyr", #data cleaning
                    #"bipartite",
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

#2 LOAD DATA--------------------------------------------------------------------

  ##2.1 LOAD GENETIC DATA----
    
    gen <- readr::read_csv(here::here("data","metabarcoding_results.csv"))
    
    #transform into a tibble to faciliate using packages like dyplr 
    gen <- as_tibble(gen)
    
  ##2.2 LOAD SEED DATA----
  
    sem <- readr::read_csv(here::here("data","seeds_28MAR2025_2.csv"))
    
    #transform into a tibble to facilitate using packages like dyplr 
    sem <- as_tibble(sem)
    
  ##2.2 LOAD SEED MORPHOSPECIES DATA----
    
    morpho <- readr::read_csv(here::here("data","morpho_4APRIL2025.csv"))
    
    #transform into a tibble to faciliate using packages like dyplr 
    morpho <- as_tibble(morpho)
    
    
#3 DATA CLEANING AND PREP-------------------------------------------------------

  ##3.1 GENETIC DATA----
    
    # cut excess rows and columns 
    
      #cut excess rows (if need be)
      gen <- gen[-(561:1246),]
      
      #cut excess columns (if need be)
      gen <- gen[,-(16:32)]
      
      #filter to only include pools that include samples with seeds, and corresponding blanks
      class(gen$pool)
      #
      gen_sem <- gen %>%
      filter(pool %in% c("P_00", "P_01", "P_02", "P_03","P_04", "P_05", "P_06",                        "P_07", "P_08", "P_09", "P_10", "P_11", "P_12"))
      
      gen_sem <- gen_sem %>%
        filter(seeds == 1 | grepl("blank", sample))
      
      #confirm number of unique samples that were sent for sequencing 
      #unique_samples <- gen %>% 
        #reframe(unique(sample)) #68 samples sent for sequencing
                                 #2 test samples
                                 #30 the first round
                                 #36 the second round 
                                
      #check number of samples successfully sequences (resulted in OTUs)
      #unique_samples_successful <- gen %>% 
        #filter(final_fate ==1) %>%
        #reframe(unique(amostra)) #59 samples successfully sequenced and 
                                  #not cut due to insufficient reads
      
      

  ##3.2 SEED & SEED MORPHOLOGY DATA----
    
    ###3.2.1 SEED DATA-----
      # cut excess rows and columns 
      
        #cut excess rows
        #sem <- sem[-(517:1023),]
        
        #cut excess columns
        sem <- sem[,-(12:19)]      
    
      #filter to only include samples with seeds
      #sem <- filter(sem, seeds == 1)
      
      # subset columns that we are interested in 
      sem_1 <- sem %>% 
        select(3, 4, 9, 10, 11)
      
      #create other subset of columns to facilitate later analyses
      #sem_2 <- sem %>% 
      #  select(3, 4, 9, 10, 11)
      
    ###3.2.2 SEED MORPHOLOGY DATA-----
     
      # cut excess rows and columns 
      
        #cut excess rows
       # morpho <- morpho[-(42:987),]
        
        #cut excess columns
        #morpho <- morpho[,-(7:9)]      
      
      # subset columns that we are interested in 
      #morpho <- morpho %>% 
        #select(2:7)
      
    ##3.2.3 MERGE SEED AND SEED MORPHOLOGY DATA----
      
      sem_morph <- merge(x=sem_1, y=morpho, by="morpho_spp", all.x =T)
      
      #sem_morph_2 <- merge(x=sem_2, y=morpho, by="morpho_spp", all.x =T)
      
      
      
      
      
      write.csv(sem_morph, here::here( "figures&results",
                                           "seed_morphospp.csv"))
      

  
  ##3.3 FILTER GENETIC DATA TO ONLY INCLUDE SEEDS WITH SAMPLES-------
    
    #only include sequenced samples that contained seeds by cross-referencing 
    #with seed list 
    #gen_sem <- gen %>%
      #filter(seeds == 1)
    
    #create a species list of all species sequenced
    spp_gen <- gen_sem %>%
      filter(seeds == 1, !is.na(seeds)) %>%
      group_by(bird_species) %>%
      summarize(n = n_distinct(sample))
    
    #write a csv for a table in manuscript
    write.csv(spp_gen, here::here( "figures&results",
                                    "spp_sequenced.csv"))
  
  ##3.4 FILTER SEED DATA TO ONLY INCLUDE SAMPLES SUCCESSFULLY SEQUENCED------ 
    
    #filter OTUs that had a final fate of 1 (were successfully sequenced)
    # 0 means it didn't make the read threshold (<0.05% of reads)
    # 2 means it didn't make the match threshold (<98% match with an organism)
    # 3 means it wasn't successfully sequenced (didn't return any OTUs)
    # 4 means it was excluded due to suspected contamination 
    
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
      group_by(species) %>%
      summarize(n = n_distinct(sample))
    
    write.csv(spp_gen_suc, here::here( "figures&results",
                                   "spp_sequenced_suc.csv"))
    
  ##3.5 MERGE MORPHOSPP DATA WITH NUMBER OF SEEDS TO CREATE SUMMARY-------
     
     #transform column no_seeds (number of seeds) into numeric
     #sem_morph_2$no_seeds <- as.numeric(sem_morph_2$no_seeds)
     
     #check that the column is now numeric
     #class(sem_morph_2$no_seeds)
     
     #create a summary table with the number of bird species in the feces of 
     #which each seed morphospp was found, as well as the number of samples, 
     #and the number of seeds of each morphospp
     
     
    
      sem_summary <- sem_morph_seq %>% 
        group_by(unique) %>%
        summarise(n_spp = n_distinct(species), 
                  n_samples = n_distinct(sample),
                  n_seed = sum(no_seeds),
                  family = first(family), genus = first(genus),
                  genus_morpho = first(genus_morpho))
    
    write.csv(sem_summary, here::here( "figures&results",
                                        "seed_morpho_summary.csv"))
    
    sem_summary_spp <- sem_morph_seq %>% 
      group_by(species) %>%
      summarise(n_morph = n_distinct(unique), 
                n_gen = n_distinct(genus))
    
    mean(sem_summary_spp$n_morph)
    
    sd(sem_summary_spp$n_morph)
    
    mean(sem_summary_spp$n_gen)
    
    sd(sem_summary_spp$n_gen)
    
    
    sem_summary_ind <- sem_morph_seq %>% 
      group_by(sample) %>%
      summarise(n_morph = n_distinct(unique), 
                n_gen = n_distinct(genus))
    
    mean(sem_summary_ind$n_morph)
    
    sd(sem_summary_ind$n_morph)
    
    mean(sem_summary_ind$n_gen)
    
    sd(sem_summary_ind$n_gen)
    
    sem_summary_plants <- sem_morph_seq %>% 
      group_by(genus) %>%
      summarise(n_spp = n_distinct(species), 
                n_ind = n_distinct(sample))
    
    mean(sem_summary_plants$n_spp)
    
    sd(sem_summary_plants$n_spp)
    
    mean(sem_summary_plants$n_ind)
    
    sd(sem_summary_plants$n_ind)
    
    #3.6 CREATE SUMMARY OF GENERA IDENTIFIED BY GENETIC DATA 
    
    #create a table of plant genera detected by metabarcoding and the number of
    #species and individuals in which it was detected in the diet
    gen_summary <- gen_sem %>%
      filter(!grepl("branco", amostra), !is.na(OTU), 
             !grepl(0, final_fate), !grepl(2, final_fate), 
             !grepl(3, final_fate), !grepl(4, final_fate)) %>%
      group_by(genus) %>%
      summarise(n_spp = n_distinct(especie), 
                n_samples = n_distinct(amostra),
                family = first(family))
    
    write.csv(gen_summary, here::here( "figures&results",
                                       "genera_gen_summary.csv"))
    
    sem_gen <- sem_morph_2 %>%
      filter(amostra %in% fil_gen_sem$amostra | grepl("branco", amostra))
    
    sum(sem_gen$no_seeds)
    
    n_distinct(sem_gen$unique)
    n_distinct(sem_gen$family)
    
# 4. SUMMARY STATISTICS AND ANALYSES---------------------------------------------
      
  ##4.1 Number of samples (with seeds) analyzed------
    
      #number of samples with seeds sent for sequencing 
      n_analyzed <- gen_sem %>%
        filter(!grepl("blank", sample)) %>%
        summarize(n_distinct(sample)) 
    
      print(n_analyzed) # 64 samples
      
      
  ##4.2 Summary of the fate of each sample----
    
      # create new column for samples successfully sequenced and used
      # the "final fate" of at least one OTU will be 1, 
      # 0 means it didn't make the read threshold,
      # 2 means it didn't make the match threshold
      # 3 means it wasn't successfully sequenced 
      # 4 means it was excluded due to suspected contamination 
      
      #create a new column for sample success (at least one OTU has the 
      # final_fate of 1)
      gen_sem <- gen_sem %>%
        group_by(amostra) %>%
        mutate(sample_success = ifelse(any(final_fate == 1), 1, 0)) %>%
        ungroup()  # Ungroup to return to the original data frame structure
      
      # now we check the fate of those that weren't successful
      sample_fate_summary <- gen_sem %>%
        group_by(amostra) %>%
        summarise(
          final_fate = ifelse(any(sample_success == 1), "successful",
                          ifelse(any(final_fate == 0), "insufficient_reads",
                            ifelse(any(final_fate == 2), "insufficient_match",
                              ifelse(any(final_fate == 3),"unsuccessful_seq", 
                                ifelse(any(final_fate == 4), "contamination"))))
                          ))
      
      #here we summarize, samples only fell into three categories:
      #   successful, unsuccessfully sequenced, or had insufficient reads
      #       while individual OTUs were deleted due to insufficient match or 
      #       contamination, those two were not motives for excluding samples
      final_fate_counts <- sample_fate_summary %>%
        filter(!grepl("branco", amostra)) %>%
        group_by(final_fate) %>%
        summarise(total_samples = n())


  ##4.2 OTU stats------
      
    # Combined read counts of filtered and unfiltered data
    # total number of paired reads
      
      # we are not going to include the blanks here, and it will be the 
      # combined read count of 12 pools: 
                # P00 (2 samples)
                # P01, 02, 03, 04, 05, 07, 08 ,09, 10, 11, 12 (6 samples each)
                # P06 and P13 are the blank pools
      
      # Quality filtering was done using a set of rules:
      
      #     1. The OTU must have > 0.05% of the total reads of the pool, if not,
      #        I assigned the final fate of 0 
      #     2. The OTU must have > 98% match with any organism (fate 2)
      #     3. If the same OTU was detected in the blank (control samples), and 
      #        it has a number of reads equal to or less than that of the 
      #        control, it is excluded (fate 4)
      
    ###4.2.1 READS AND OTUs BEFORE QUALITY FILTERING-----
      
        # total number of reads per OTU, filtering out the blank samples and 
        # the OTUs that have a value NA because they weren't successfully 
        # sequenced
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
        
        #the number of OTUs produced
        n_OTUs <- t_reads_OTU %>%
            summarise(n_distinct(OTU))
        
        print(n_OTUs)
      
      ###4.2.2 READS AND OTUs AFTER QUALITY FILTERING-----
        
        fil_gen_sem <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(OTU) %>%
          summarise(t_reads = sum(no_reads_OTU))
        # total number of reads per OTU, filtering out the blank samples and 
        # the OTUs that have a value NA because they weren't successfully 
        # sequenced and now, the OTUs that were excluded due to insufficient
        # reads or match 
        # t_reads_OTU_fil <- gen_sem %>%
        #   filter(!grepl("branco", amostra), !is.na(OTU), 
        #          !grepl(0, final_fate), !grepl(2, final_fate), 
        #          !grepl(3, final_fate), !grepl(4, final_fate)) %>%
        #   group_by(OTU) %>%
        #   summarise(t_reads = sum(no_reads_OTU))
        
        genera <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          summarize(n = n_distinct(plant_genus))
        
        print(genera)
        
        #the sum of filtered reads
        sum(fil_gen_sem$t_reads)
        
        #the mean number of reads per OTU
        mean(fil_gen_sem$t_reads)
        
        #the standard deviation of reads
        sd(fil_gen_sem$t_reads)
        
        #the number of OTUs produced
        n_OTUs_fil <- fil_gen_sem %>%
          summarise(n_distinct(OTU))
        
        print(n_OTUs_fil)
        
        #number of reads per spp (not adjusted for sample size)
        t_reads_fil_spp <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(bird_species) %>%
          summarise(n_reads_spp = sum(no_reads_OTU))
        
        mean(t_reads_fil_spp$n_reads_spp)
        
        sd(t_reads_fil_spp$n_reads_spp)
        
        #adjust for sample size 
        colnames(spp_gen_suc)[1] <- "bird_species"
        spp_gen_suc
        
        reads_fil_spp_ind <- merge(x=t_reads_fil_spp, y=spp_gen_suc, 
                           by="bird_species", all.x =T)
        
        reads_fil_spp_ind$n_reads_ind <- (reads_fil_spp_ind$n_reads_spp / 
                                         reads_fil_spp_ind$n)
        
        mean(reads_fil_spp_ind$n_reads_ind)
        
        sd(reads_fil_spp_ind$n_reads_ind)
        
        #number of reads per ind 
        t_reads_fil_ind <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(sample) %>%
          summarise(n_reads_ind = sum(no_reads_OTU))
        
        mean(t_reads_fil_ind$n_reads_ind)
        
        sd(t_reads_fil_ind$n_reads_ind)
        
        
        #number of OTU per spp (not adjusted for sample size)
        t_reads_OTU_fil_spp <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(bird_species) %>%
          summarise(n_OTU_spp = n_distinct(OTU))
        
        mean(t_reads_OTU_fil_spp$n_OTU_spp)
        
        sd(t_reads_OTU_fil_spp$n_OTU_spp)
        
        #number of OTU per ind (not adjusted for sample size)
        t_reads_OTU_fil_ind <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(sample) %>%
          summarise(n_OTU_ind = n_distinct(OTU))
        
        mean(t_reads_OTU_fil_ind$n_OTU_ind)
        
        sd(t_reads_OTU_fil_ind$n_OTU_ind)
        
        #number of genera per spp
        t_reads_gen_fil_spp <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(bird_species) %>%
          summarise(n_gen_spp = n_distinct(plant_genus))
        
        mean(t_reads_gen_fil_spp$n_gen_spp)
        
        sd(t_reads_gen_fil_spp$n_gen_spp)
        
        #number of genera per ind
        t_reads_gen_fil_ind <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(sample) %>%
          summarise(n_gen_ind = n_distinct(plant_genus))
        
        mean(t_reads_gen_fil_ind$n_gen_ind)
        
        sd(t_reads_gen_fil_ind$n_gen_ind)
        
      ###4.2.3 OTU AND READ FATE SUMMARY------
        
        OTU_fate_summary <- gen_sem %>%
          filter(!grepl("branco", amostra),!is.na(OTU)) %>%
          group_by(OTU) %>%
          summarise(
            final_fate = ifelse(any(final_fate == 1), "successful",
                          ifelse(any(final_fate == 0), "insufficient_reads",
                            ifelse(any(final_fate == 2), "insufficient_match",
                              ifelse(any(final_fate == 3),"unsuccessful_seq", 
                                ifelse(any(final_fate == 4), "contamination"))))
            ), no_reads_OTU = no_reads_OTU)
        
        OTU_fate_counts <- OTU_fate_summary %>%
          group_by(final_fate) %>%
          summarise(total_OTUs = n_distinct(OTU), n_reads = sum(no_reads_OTU))
        
  ## 4.3 DETECTIONS: PLANT FAMILIES AND GENERA--------
        
        n_plant_fam <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(plant_family) %>%
          summarise(n_plant_fam = n())
        
        n_distinct(n_plant_fam$plant_family)
        
        n_plant_gen <- gen_sem %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(plant_genus) %>%
          summarise(n_spp = n_distinct(bird_species), 
                    n_ind = n_distinct(sample))
        
        mean(n_plant_gen$n_spp)
        
        sd(n_plant_gen$n_spp)
        
        mean(n_plant_gen$n_ind)
        
        sd(n_plant_gen$n_ind)
        
#       3j. Number of detections / OTU --> accumulation curve----
        
# 5 TRIPARTITE NETWORK ANALYSIS-------------------------------------------------


  #5.1 MAKE AN INTERACTION MATRIX WITH GENETIC DATA----
          
        #make a matrix
        gen_mat <- gen_sem_seq %>%
          filter(!grepl("blank", sample), !is.na(OTU)) %>%
          filter(final_fate == 1)  %>%
          group_by(sample, plant_genus) %>%  # Group by samples and genera
          summarise(occurrences = n(), .groups = 'drop') %>%  # Count occurrences for each sample-genus pair
          pivot_wider(names_from = plant_genus, values_from = occurrences, values_fill = 0) %>%  # Wide format: genera as columns
          column_to_rownames(var = "sample")
        
        #write.csv(gen_mat, here::here( "figures&results",
                                        #"gen_mat.csv"))
        gen_mat[gen_mat > 1] <- 1
        
        #first, we need to make two matrices, one for the genetic data and
        #one for the morphological data, but they need to be the same 
        #dimensions, so we are first going to make the row names into the first
        #column so that we can use the melt function to transform it into a list 
        
        gen_mat_with_rowname <- rownames_to_column(gen_mat, var = "Sample")
        
        #transform the genetic data into a list 
        gen_list <- melt(gen_mat_with_rowname)
        
        #rename the columns, in this case, "Weight" is the number of times that 
        #the genus was detected in that sample
        colnames(gen_list) <- c("Sample", "Genus", "Weight")
        
        #but! we want presence absence data, because sometimes multiple OTUs 
        #that we classified as the same genus were present in the sample, 
        #making it look like, for example, we detected Miconia 4 times in the 
        #same sample, so we are going to transform any number >1 into 1
        gen_list$Weight[gen_list$Weight > 1] <- 1
        
        df_gen <- subset(gen_list, Weight >= 1)

        desired_order_samples <- c("040", "050", "052", "084",
                                   "070","043", "053", "272", 
                                   "432", "073","419", 
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
                                   "290", 
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
                                "Myrciaria","Myrcia", "Psidium",
                                "Morus", "Rubus",  "Siparuna"
        )
        
        
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
              "Syagrus" = "#000000",
              "Myrciaria" = "#000000"# Black
              
            ))

        ggsave(
          filename = "gen_graph.svg", # File name
          plot = last_plot(),             # Optional: explicitly specify the plot
          width = 4,                      # Width in inches
          height = 16,                     # Height in inches
          dpi = 720                       # Resolution in dots per inch
        )
        
        write.csv(df_gen, here::here( "figures&results",
                                       "df_gen.csv"))
        
  #5.2 MAKE AN INTERACTION MATRIX WITH SEED DATA----
        
      
         #by genus, with genera as columns and samples as rows 
         seed_mat <- sem_morph_seq %>%
           group_by(sample, genus) %>%  # Group by samples and genera
           summarise(occurrences = n(), .groups = 'drop') %>%  # Count occurrences for each sample-genus pair
           pivot_wider(names_from = genus, values_from = occurrences,
                       values_fill = 0) %>%  # Wide format: genera as columns
           column_to_rownames(var = "sample")
        
        #now we are going to do the same thing with the morphology data 
        seed_mat_with_rowname <- rownames_to_column(seed_mat, 
                                                    var = "Sample")
        seed_list <- melt(seed_mat_with_rowname)
        
        colnames(seed_list) <- c("Sample", "Genus", "Weight")
        
        seed_list$Weight[seed_list$Weight > 1] <- 1
        
        
        write.csv(seed_mat, here::here( "figures&results",
                                         "seed_mat.csv"))
        
        df_seed <- subset(seed_list, Weight >= 1)
        
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
                                 "Lantana",
                                 "Solanum", "Miconia",
                                 
                                 
                                 "Allophylus", "Maclura", "Myrtaceae_2",
                                 
                                 
                                 "Stylogyne", 
                                 "Fabaceae_1",
                                 "Myrcia", "unknown_1","Psidium","Morus",
                                 
                                 
                                 
                                 
                                 "Rubus","Struthanthus","Pera",
                                 "Siparuna"
        )
        

        df_seed$Sample <- factor(df_seed$Sample, levels = desired_order_samples)
        
        df_seed$Genus <- factor(df_seed$Genus, levels = desired_order_seed)
        
        
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
        
        ggsave(
          filename = "morph_graph.svg", # File name
          plot = last_plot(),             # Optional: explicitly specify the plot
          width = 4,                      # Width in inches
          height = 16,                     # Height in inches
          dpi = 720                       # Resolution in dots per inch
        )
        
        write.csv(df_seed, here::here( "figures&results",
                                        "df_seed.csv"))
#6 JACCARD INDEX AND RECALL-----------------------------------------------------
        
  #6.1 CREATE MATRICIES OF THE TWO DATA SETS WITH THE SAME DIMENSIONS----
        
      #6.1.1 GENETIC DATA----
        
        #first, we need to make two matrices, one for the genetic data and
        #one for the morphological data, but they need to be the same 
        #dimensions, so we are first going to make the row names into the first
        #column so that we can use the melt function to transform it into a list 
        
        gen_mat_with_rowname <- rownames_to_column(gen_mat, var = "Sample")
        
        #transform the genetic data into a list 
        gen_list <- melt(gen_mat_with_rowname)
        
        #rename the columns, in this case, "Weight" is the number of times that 
        #the genus was detected in that sample
        colnames(gen_list) <- c("Sample", "Genus", "Weight")
        
        #but! we want presence absence data, because sometimes multiple OTUs 
        #that we classified as the same genus were present in the sample, 
        #making it look like, for example, we detected Miconia 4 times in the 
        #same sample, so we are going to transform any number >1 into 1
        gen_list$Weight[gen_list$Weight > 1] <- 1
        
      #6.1.2 SEED MORPHOLOGY DATA----
        
        #now we are going to do the same thing with the morphology data 
        seed_mat_with_rowname <- rownames_to_column(seed_mat, 
                                                    var = "Sample")
        seed_list <- melt(seed_mat_with_rowname)
        
        colnames(seed_list) <- c("Sample", "Genus", "Weight")
        
        seed_list$Weight[seed_list$Weight > 1] <- 1
        
        gen_list_summary <- gen_list %>%
          group_by(Genus) %>%
          summarise(n = sum(Weight))
        
        ordered <-  gen_list_summary[order(-gen_list_summary$n), ]
        
        print(ordered, n = 41)
        
      #6.1.3 FIND THE MISSING DATA TO MAKE THE MATRICES THE SAME DIMENSIONS-----
        
        #6.1.3.1 FIND DATA MISSING FROM MORPHOLOGY DATA SET-----
        
        # Find the samples missing from the seed morphology data set that are 
        #in the genetic data set
        missing_from_morph <- gen_list %>% 
          anti_join(seed_list, by = "Genus")
        
        #make their presence 0 (they were not present, but we need the names
        #so the matrices are the same dimensions)
        missing_from_morph$Weight <- 0
        
        # Combine both data frames, including missing data
        seed_complete <- bind_rows(seed_list, missing_from_morph)
        
        #6.1.3.2 FIND DATA MISSING FROM GENETIC DATA SET-----
        
        # Find missing samples in gen that are in seeds
        missing_from_gen <- seed_list %>% 
          anti_join(gen_list, by = "Genus")
        
        #make their presence 0 (they were not present, but we need the names
        #so the matrices are the same dimensions)
        missing_from_gen$Weight <- 0
        
        # Combine both data frames, including missing data
        gen_complete <- bind_rows(gen_list, missing_from_gen)
        
        #6.1.3.3 RETURN THE LISTS TO MATRICES-----
        
        
        #a complete matrix of the genetic data, with the same number of rows 
        #and columns as the morphological data now that we included the missing
        #data 
        genetic_matrix <- gen_complete %>%
          pivot_wider(names_from = Genus, values_from = Weight) %>%
          column_to_rownames(var = "Sample")
        
        #a complete matrix of the morphological data, with the same number of 
        #rows and columns as the genetic data now that we included the missing
        #data 
        morpho_matrix <- seed_complete %>%
          pivot_wider(names_from = Genus, values_from = Weight) %>%
          column_to_rownames(var = "Sample")

        #put the matrices in the same order
        genetic_matrix <- genetic_matrix[, order(colnames(genetic_matrix))]
        
        morpho_matrix <- morpho_matrix[, order(colnames(morpho_matrix))]
        
        #make the row names the first column so we can run the Jaccard Index 
        #and recall 
        genetic_matrix_w_rowname <- rownames_to_column(genetic_matrix, 
                                                    var = "Sample")
        
        morpho_matrix_w_rowname <- rownames_to_column(morpho_matrix, 
                                                       var = "Sample")
      
        
  #6.2 JACCARD--------------
        
        # Function to calculate Jaccard Index for two binary vectors
        calculate_jaccard <- function(x, y) {
          intersection <- sum(x == 1 & y == 1)
          union <- sum(x == 1 | y == 1)
          return(intersection / union)
        }
        
        # Apply the Jaccard Index calculation for each sample
        overlap <- sapply(1:nrow(genetic_matrix_w_rowname), function(i) {
          calculate_jaccard(genetic_matrix_w_rowname[i, -1], 
                            morpho_matrix_w_rowname[i, -1])
        })
        
        # Create a data frame with the overlap results
        result_jaccard <- data.frame(
          Sample = genetic_matrix_w_rowname$Sample,
          Jaccard_Index = overlap
        )
        
        # View the result
        print(result_jaccard)
        
        mean(result_jaccard$Jaccard_Index)
        
        sd(result_jaccard$Jaccard_Index)
        
  #6.3 RECALL (morphology as true data) -------------------
        # Combine the lists into one data frame
        combined_list <- bind_rows(df_seed, df_gen)
        
        # Group by Sample_Number and count unique genera per sample
        genera_per_sample <- combined_list %>%
          group_by(Sample) %>%
          summarise(Total_Unique_Genera = n_distinct(Genus))
        
        #how well the genetic method captures the genera identified by morphology 

        calculate_recall <- function(morph_data, genetic_data) {
          true_positives <- sum(morph_data == 1 & genetic_data == 1)
          false_negatives <- sum(morph_data == 1 & genetic_data == 0)
          
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        }
        
        # Apply recall calculation for each sample
        recall_values <- sapply(1:nrow(genetic_matrix_w_rowname), 
                                function(i) {
          calculate_recall(morpho_matrix_w_rowname[i, -1], 
                           genetic_matrix_w_rowname[i, -1])
        })
        
        # Create a data frame with the recall results
        recall_result_morph <- data.frame(
          Sample = genetic_matrix_w_rowname$Sample,
          Recall = recall_values
        )
        
        mean(recall_result_morph$Recall)
        
        
  #6.4 RECALL (genetic as true data) -------------------
        
        
        
  #RECALL AGAIN
        
        df_gen$Genus <- as.character(df_gen$Genus)
        df_list$Genus <- as.character(df_list$Genus)
        
        merged <- merge(df_gen, df_list, by = "Sample", 
                             suffixes = c("_gen", "_seed"))
        
        
        
        # Add a column to check if genera are the same
        merged$Genus_Match <- merged$Genus_gen == merged$Genus_seed
        
        merged_2 <- merged %>%
          group_by(Sample) %>%
          summarise(n_match = sum(Genus_Match))
        
        
        
        genera_summary_1 <- 
          merge(genera_summary, merged_2, by = "Sample", all = FALSE)
        
        print(mean(genera_summary_1$r_gen))
        
        
        print(mean(genera_summary$r_morpho))
        
        class(genera_summary$n_match)
        
        genera_summary_1$n_match_gen <- (genera_summary_1$n_match)/
                                        (genera_summary_1$Morphological_Count)
        
        genera_summary_1$n_match_morph <- (genera_summary_1$n_match)/
          (genera_summary_1$Genetic_Count)
        
        genera_summary_1$p_method_morph <- (genera_summary_1$Morphological_Count)/
          (genera_summary_1$Total_Unique_Genera)
        
        mean(genera_summary_1$p_method_gen)
        
        genera_summary_1$p_method_gen <- (genera_summary_1$Genetic_Count)/
          (genera_summary_1$Total_Unique_Genera)
        
        genera_summary_1$recall_dif <- (genera_summary_1$n_match_morph - 
                                        genera_summary_1$n_match_gen)
        
        mean(genera_summary_1$p_method_morph)
        
        print(mean(genera_summary_1$n_match_gen))
        
        print(mean(genera_summary_1$n_match_morph))
        
        print(sum(genera_summary_1$n_match_gen))
        
        print(sum(genera_summary_1$n_match_morph))
        
        print(sum(genera_summary_1$eval_method))
        
        print(sum(genera_summary_1$Morphological_Count))
        
        print(sum(genera_summary_1$Genetic_Count))
        
      genera_summary_1$Comparison <- ifelse(genera_summary_1$Genetic_Count 
                                            > genera_summary_1$Morphological_Count, "Greater",
                                ifelse(genera_summary_1$Genetic_Count 
                                       < genera_summary_1$Morphological_Count, "Less", "Equal"))
      comparison_summary <- genera_summary_1 %>%
        group_by(Comparison) %>%
        summarise(
          Count = n(),
          Equal_but_different = sum(Comparison == "Equal" & n_match < Total_Unique_Genera)
        )
     #how well the morphological method captures the genera identified by 
     #genetics 
        
    # Define the recall calculation function
        calculate_recall <- function(genetic_data, morph_data) {
          true_positives <- sum(genetic_data == 1 & morph_data == 1)
          false_negatives <- sum(genetic_data == 1 & morph_data == 0)
          
          # Handle cases where denominator is 0
          if ((true_positives + false_negatives) == 0) {
            return(NA)  # Or return 0, depending on how you want to handle it
          } else {
            recall <- true_positives / (true_positives + false_negatives)
            return(recall)
          }
        }
        
        # Apply recall calculation for each sample
        recall_values <- sapply(1:nrow(morpho_matrix_w_rowname), 
                                function(i) {
                                  calculate_recall(genetic_matrix_w_rowname[i, -1], 
                                                   morpho_matrix_w_rowname[i, -1])
                                })
        
        # Create a data frame with the recall results
        recall_result_gen <- data.frame(
          Sample = morpho_matrix_w_rowname$Sample,
          Recall = recall_values
        )
        
        # Calculate the mean of recall values
        mean_recall <- mean(recall_result_gen$Recall)
        
        # Print the mean recall
        print(mean_recall)
        
    #6.5 RECALL PAIRED T-TEST------
        
        s_t_morph <- shapiro.test(genera_summary_1$p_method_morph)
        
        s_t_gen <- shapiro.test(genera_summary_1$p_method_gen)
        
        print(s_t_morph)
        
        print(s_t_gen)
        
        recall_result_gen$log <- log(genera_summary_1$p_method_morph + 0.001)
        
        recall_result_morph$log <- log(genera_summary_1$p_method_gen + 0.001)
        
        s_t_morph_log <- shapiro.test(recall_result_gen$log)
        
        print(s_t_morph_log)
        
        t_test_result <- t.test(genera_summary_1$p_method_morph, 
                                genera_summary_1$p_method_gen, 
                         paired = TRUE)
        
        print(t_test_result)
        
        wilcox_test_result <- wilcox.test(genera_summary_1$p_method_morph, 
                                          genera_summary_1$p_method_gen, 
                                          paired = TRUE)
        
        print(wilcox_test_result)
        
        
        recall_combined <- data.frame(
          Recall = c(recall_result_gen$Recall, recall_result_morph$Recall),
          Method = factor(c(rep("Genetic", length(recall_result_gen$Recall)), 
                            rep("Morphological", length(recall_result_morph$Recall))))
        )
        
        library(coin)
        # Perform a permutation test
        wilcox_test(Recall ~ Method, data = recall_combined, 
                    distribution = "exact")
        
        
        
    #6.6 RECALL DIFFERENCE-------------
        
        recall_difference <- (recall_result_morph$Recall - recall_result_gen$Recall)
        
        #sample_simple <- (1:54)
        
        recall_difference_df <- data.frame(recall_result_morph$Sample, 
                                           recall_difference)
        
        colnames(recall_difference_df) <- c("Sample", "recall_difference")
        

        #6.6.1 PLOTTING RECALL DIFFERENCE VALUES AS COLOR GRADIENT ------
        
        color_gradient_recall <- scale_color_gradient2(
          low = "#4100ff",    # Darkest blue at -1
          mid = "#a240bf",  # Purple at 0
          high = "#ff0045",    # Darkest red at 1
          midpoint = 0,    # The center point of the gradient
          limits = c(-1, 1) # The range of your values
        )
        
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
        
        # Reorder the dataframe
       genera_summary_1$Sample <- 
          factor(genera_summary_1$Sample, levels = order_samples) # Apply the order
        
        recall_difference_df_ordered <- genera_summary_1[order(genera_summary_1$Sample), ] # Reorder the rows
        
        
        
        
        #visualization
        
        ggplot(recall_difference_df_ordered, aes(x = eval_method, y = Sample)) +
          geom_segment(aes(x = 0, xend = eval_method, y = Sample, yend = Sample), color = "grey") +  # Horizontal line
          geom_point(aes(x = eval_method, y = Sample), 
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
        
        
        #gray
        ggplot(recall_difference_df_ordered, aes(x = recall_difference, y = Sample)) +
          geom_segment(aes(x = 0, xend = recall_difference, y = Sample, yend = Sample), color = "black") +  # Horizontal line
          geom_point(aes(x = recall_difference, y = Sample, color = recall_difference), 
                     size = 3) +  # Dot at the end
          geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Vertical line at 0
          labs(title = "Difference Between Genetic and Morphological Methods for Each Sample",
               x = "Difference (Genetic - Morphological)",
               y = "Sample") +
          scale_color_gradient(low = "black", high = "lightgray") +  # Color gradient between black and white
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
        
        # Plot
        # ggplot(recall_difference_df, aes(x = 0, y = sample_simple, 
        #                                  color = recall_difference)) +
        #   geom_point(size = 5) +
        #   color_gradient_recall +  # Add the gradient
        #   theme_minimal() +
        #   theme(
        #   panel.grid = element_blank(),  # Remove grid lines
        #   axis.title.x = element_blank(),  # Remove x-axis title
        #   axis.text.x = element_blank(),   # Remove x-axis labels
        #   axis.ticks.x = element_blank()   # Remove x-axis ticks
        # ) +
        #   scale_y_continuous(name = "Sample", 
        #                      expand = expansion(mult = c(0.001, 0.001)))
        
        
    #6.7 MARINA RECALL------
        
        # Assuming genetic_matrix_w_rowname and morpho_matrix_w_rowname have 0/1 values
        # and the first column is the sample names
        
        # Function to calculate the number of genera identified by each method
        calculate_genera_counts <- function(genetic_data, morph_data) {
          
          # Create a data frame to store the results
          genera_df <- data.frame(
            Sample = genetic_data[, 1],  # Assuming first column is Sample names
            Genetic_Count = rowSums(genetic_data[, -1]),   # Sum of genera identified by genetic method
            Morphological_Count = rowSums(morph_data[, -1])  # Sum of genera identified by morphological method
          )
          
          # Calculate the total number of unique genera (union of genetic and morphological)
          genera_df$Total_Unique_Genera <- sapply(1:nrow(genetic_data), function(i) {
            sum((genetic_data[i, -1] + morph_data[i, -1]) > 0)  # Union of both methods for each sample
          })
          
          return(genera_df)
        }
        
        # Apply the function to your matrices
        genera_summary <- calculate_genera_counts(genetic_matrix_w_rowname, morpho_matrix_w_rowname)
        
        # View the result
        print(genera_summary)
        
        #calculate recall using morphological ids as "true" dataset
        genera_summary$r_morpho <- genera_summary$Genetic_Count/
                                    genera_summary$Total_Unique_Genera
        
        #calcuate recall using morphological ids as "true" dataset
        genera_summary$r_gen <- genera_summary$Morphological_Count/
                                  genera_summary$Total_Unique_Genera
        
        #now I'm going to subtract r_gen from r_morph, if the value is 0,
        #then they identified the same stuff, if the value is greater than 
        #1, then the genetic method identified more, and if it is 
        genera_summary$eval_method <- genera_summary$r_morpho - 
                                               genera_summary$r_gen
        
        #visualization 
        
        ggplot(genera_summary_1, aes(x = recall_dif)) +
          geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +  # Adjust binwidth as needed
          geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a vertical line at 0
          labs(title = "Difference Between Genetic and Morphological Methods",
               x = "Difference (Genetic - Morphological)",
               y = "Number of Samples") +
          theme_minimal()
        
        # Create a lollipop chart
        ggplot(genera_summary_1, aes(x = recall_dif, y = Sample)) +
          geom_segment(aes(x = 0, xend = eval_method, y = Sample, yend = Sample), color = "grey") +  # Horizontal line
          geom_point(aes(x = eval_method, y = Sample), color = "blue", size = 3) +  # Dot at the end
          geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Vertical line at 0
          labs(title = "Difference Between Genetic and Morphological Methods for Each Sample",
               x = "Difference (Genetic - Morphological)",
               y = "Sample") +
          theme_minimal() +
          theme(panel.grid = element_blank(),  # Remove all gridlines
                axis.text.y = element_text(size = 8)) +
          xlim(-3, 3) # Adjust size if needed  # Adjust size if needed
#7 VENN DIAGRAM-----------------------------------------------------------------
        
        class(gen_list$Weight)
        
        list_detections <- gen_list %>%
          filter(Weight == 1) %>%
          distinct(Genus, .keep_all = TRUE)
        
        list_det_onlygen <- list_detections[-c(34,11), ]

        
        list_seeds <- seed_list %>%
          filter(Weight == 1) %>%
          distinct(Genus, .keep_all = TRUE)
        
        n_distinct(list_seed_onlygen$Genus)
        
        writeClipboard(as.character(list_det_onlygen$Genus))
        list_seed_onlygen <- list_seeds[-c(26:29, 22, 20, 18), ]

        venn.plot <- venn.diagram(
          x = list(Genetic = list_det_onlygen$Genus, 
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
        
        # Display the Venn diagram
        grid.draw(venn.plot)       
        
####------ADJUSTED RECALL --> fuzzy values--------------------------------------
        
        
        gen_mat_fuzzy <- gen_sem %>%
          filter(!grepl("branco", amostra), 
                 !grepl(0, final_fate), 
                 !grepl(4, final_fate)) %>%
          group_by(amostra, genus) %>%  # Group by samples and genera
          summarise(occurrences = n(), .groups = 'drop') %>%  # Count occurrences for each sample-genus pair
          pivot_wider(names_from = genus, values_from = occurrences, values_fill = 0) %>%  # Wide format: genera as columns
          column_to_rownames(var = "amostra")
        
        adjusted_recall <- function(genetic_pred, morph_pred, threshold = 0.3) {
          # Apply threshold to genetic data
          adjusted_genetic <- ifelse(genetic_pred >= threshold, 1, 0)
          
          true_positive <- sum(adjusted_genetic & morph_pred)  # Correctly identified by both
          total_true <- sum(morph_pred)  # Total identified by morphological method
          return(true_positive / total_true)
        }
        
        # Calculate recall with adjustment
        recall_gen_adjusted <- adjusted_recall(genetic_pred, morph_pred)
        print(paste("Adjusted Recall (Genetic captures Morphological):", recall_gen_adjusted))
        
        