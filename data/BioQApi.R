
library(plumber)

library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(rredlist)
library(vegan)
library(jsonlite)
library(FD)
library(igraph)
library(bipartite)
library(tidyr)
library(ggplot2)
library(terra)
library(raster)
library(fundiversity)
library(randomForest)
library(DT)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

#* @apiTitle Plumber Example API
#* @apiDescription Plumber example description.


#* @post /upload
#* Upload a file and process data
#* @consumes multipart/form-data
#* @param file:file The CSV file to upload
#* @return  The summary data of the uploaded CSV

function(file) {
  
  
  data_table <- as.data.table(file)
  
  
  Species_Summary <- fread(file.path(file))

  Species_Summary <- Species_Summary[,c("Sample", "Species", "Abundance")]
  Trait_list <- fread(file.path("CEFAS_Trait_DataBase.csv"))
  
  Species_Summary <- subset(Species_Summary, Abundance != 0)
  
  Species_Summary$V1<-NULL
  
  #========================================================================
  #
  # TAXONOMIC INDICI CALCULATIONS
  #
  # by Lars O. Mortensen, DHI June 2022
  #
  #========================================================================
  
  Taxonomic_indicies<-c()
  
  #=============================
  # Number of Individuals and Species
  #=============================
  
  #N
  Species.N <- ddply(Species_Summary, .(Sample), summarise, N = sum(Abundance))
  
  Species.N <- Species_Summary %>% 
    group_by(Sample) %>% 
    summarize(N = sum(Abundance))
  
  # S
  Species.S <- Species_Summary %>% 
    group_by(Sample, Species) %>% 
    summarize(S = sum(Abundance))
  
  Species.S <- Species.S %>% 
    group_by(Sample) %>% 
    summarize(S = sum(S > 0, na.rm = TRUE))
  
  #=============================================
  # Species richness (d, Margaleff) d = (S - 1) / ln N
  #=============================================
  
  # S = number of species, N = total number of individuals in the sample
  Richness <- merge(Species.S, Species.N, by = c("Sample"), all = TRUE )
  Richness$d <- (Richness$S-1) / log(Richness$N)
  
  #=============================================
  # Diversity (Simpson index (1-lambda))
  #=============================================
  
  #Making Matrix
  Div_matrix <- acast(Species_Summary, Sample ~ Species , value.var='Abundance', fun.aggregate=sum, margins=FALSE)
  
  #Calculating simpson diveristy
  Simp <- vegan::diversity(Div_matrix, "simpson")
  Simp <- data.frame(Simp)
  Simp$Sample <- rownames(Simp)
  
  #=============================================
  # Diversity (Shannon)
  #=============================================
  
  #Calculating shannon diveristy
  Shan <- vegan::diversity(Div_matrix, "shannon")
  Shan <- data.frame(Shan)
  Shan$Sample <- rownames(Shan)
  
  #=============================================
  # Aggregating last diversity table
  #=============================================
  
  Diversity_table <- merge(Richness, Simp, by = "Sample")
  Diversity_table <- merge(Diversity_table, Shan, by = "Sample")
  
  #=============================================
  # Calculating Pielous index
  #=============================================
  
  Diversity_table$Pielou <- Diversity_table$Shan/log(Diversity_table$S)
  
  #=============================================
  # Non-indigenous species
  ## @knitr NIS_overivew
  #=============================================
  
  # #Searching EASIN
  # Species_names <- unique(Species_Summary$Species)
  # 
  # Species_names <- unique(Species_Summary$Species)
  # Species_names_EOL <- Species_names
  # Species_names <- gsub(" ", "%20", Species_names)
  # 
  # Species_status <- data.frame()
  # cal <- 0
  # for(i in Species_names){
  #   AphiaID <- i
  #   url <- sprintf("https://easin.jrc.ec.europa.eu/apixg/catxg/term/%s", AphiaID);
  # 
  #   #Get the actual data from the URL
  #   classificationTree <- fromJSON(url)
  #   
  #   if(length(classificationTree) == 1){
  #     Name <- i
  #     IsEUConcern <- FALSE
  #     needed <- as.data.frame(cbind(Name, IsEUConcern))
  #     needed$Name <- gsub("%20", " ", needed$Name)
  #   }
  #   
  #   if(length(classificationTree) > 1){
  #     needed <- classificationTree[,c("Name", "IsEUConcern")]
  #     needed <- needed[grepl(i, needed$Name, fixed = TRUE),]
  #     if(nrow(needed) > 0){needed$Name <- i}
  #     if(nrow(needed) == 0){
  #       needed <- classificationTree[,c("Name", "IsEUConcern")]
  #       needed <- needed[1,]
  #       needed$Name <- i
  #       needed$IsEUConcern <- FALSE
  #       }
  #     
  #   }
  #   
  #   needed$NIS <- ifelse(needed$IsEUConcern == "FALSE", 0, 1)
  #   needed <- ddply(needed, .(Name), summarise, NIS = sum(NIS))
  #   needed$Name <- gsub("%20", " ", needed$Name)
  #   
  #   Species_status <- rbind(Species_status, needed)
  #   
  #   cal <- cal+1
  #   print(cal/length(Species_names)*100)
  # }
  # 
  # names(Species_status) <- c("Species", "NIS")
  # Alien <- merge(Species_Summary, Species_status, by = "Species", all = TRUE)
  #  
  # Alien$NIS <- Alien$NIS * Alien$Abundance
  # Alien_sub <- ddply(Alien, .(Sample), summarise, NIS = sum(NIS))
  # 
  # Diversity_table <- merge(Diversity_table, Alien_sub, by = c("Sample"), all.x = TRUE)
  # 
  # #=============================================
  # # Estimating rarity
  # #=============================================
  # 
  # data_species <- Species_Summary
  # 
  # Rarity <- data.frame()
  # cal <- 0
  # for(r in unique(Species_Summary$Species)){
  #   data <- rl_history(r, key = "3d359a007f87594b24edd35a3d8a53bfb7a7ca81657eb97ea8f4ee2d958d82dd")
  #   status <- data$result$category
  #   status <- ifelse(length(status)> 0, status, "NON")
  #   result <- cbind(r, status)
  #   Rarity <- rbind(Rarity, result)
  #   
  #   cal <- cal+1
  #   print(cal/length(unique(Species_Summary$Species))*100)
  # }
  # 
  # names(Rarity) <- c("Species", "status")
  # 
  # data_species <- merge(data_species, Rarity, by = "Species")
  # data_species$status.factor <- ifelse(data_species$status == "NON", 0, 1)
  # data_species$status.abund <- ifelse(data_species$status == "NON", 0,data_species$Abundance)
  # 
  # Rarity_species <- ddply(data_species, .(Sample), summarise, n.rare = sum(status.factor), a.rare = sum(status.abund))
  # Diversity_table <- merge(Diversity_table, Rarity_species, by = "Sample", all = TRUE)
  
  Taxonomic_indicies <- base::rbind(Taxonomic_indicies, Diversity_table)
  
  
  print("1_BioQ - Done")
  

  
  
  #========================================================================
  #
  # FUNCTIONAL INDICIE CALCULATIONS
  #
  # by Lars O. Mortensen, DHI Feb 2023
  #
  #========================================================================
  
  
  Platform_year_genus <- ddply(Species_Summary, .(Sample, Species), summarise, 
                               Abundance = sum(Abundance, na.rm = TRUE))
  
  names(Trait_list)[names(Trait_list) == '_Order'] <- 'Order'
  
  #Trait_simplification
  Trait_list_Family <- subset(Trait_list, Genus == "")
  Trait_list_Order <- subset(Trait_list, Genus == "" & Family == "")
  Trait_list_Class <- subset(Trait_list, Genus == "" & Family == "" & Order == "")
  Trait_list_Phylum <- subset(Trait_list, Genus == "" & Family == "" & Order == "" & Class == "")
  
  #Subsertting data
  
  Trait_genus <- subset(Platform_year_genus, Species %in% Trait_list$Genus)
  Trait_genus$Genus <- Trait_genus$Species
  Family_genus <- subset(Platform_year_genus, Species %in% Trait_list_Family$Family)
  Family_genus$Family <- Family_genus$Species
  Order_genus <- subset(Platform_year_genus, Species %in% Trait_list_Order$Order)
  Order_genus$Order <- Order_genus$Species
  Class_genus <- subset(Platform_year_genus, Species %in% Trait_list_Class$Class)
  Class_genus$Class <- Class_genus$Species
  Phylum_genus <- subset(Platform_year_genus, Species %in% Trait_list_Phylum$Phylum)
  Phylum_genus$Phylum <- Phylum_genus$Species
  Leftover <- subset(Platform_year_genus, !Species %in% c(Trait_list$Genus, Trait_list_Family$Family,
                                                          Trait_list_Order$Order, Trait_list_Class$Class, Trait_list_Phylum$Phylum))
  
  Numbers_iden <- as.data.frame(cbind(length(unique(Trait_genus$Genus)), 
                                      length(unique(Family_genus$Family)), 
                                      length(unique(Order_genus$Order)), 
                                      length(unique(Class_genus$Class)), 
                                      length(unique(Phylum_genus$Phylum)), 
                                      length(unique(Leftover$Genus))))
  
  
  
  names(Numbers_iden) <- c("Genus", "Familiy", "Order", "Class", "Phylum", "Unident")
  
  Trait_genus <- merge(Trait_genus, Trait_list, by = "Genus", all.x = TRUE)
  Family_genus <- merge(Family_genus, Trait_list_Family, by = "Family", all.x = TRUE)
  Family_genus$Genus <- Family_genus$Species
  Family_genus <- Family_genus[,names(Trait_genus)]
  Order_genus <- merge(Order_genus, Trait_list_Order, by = "Order", all.x = TRUE)
  Order_genus$Genus <- Order_genus$Species
  Order_genus <- Order_genus[,names(Trait_genus)]
  Class_genus <- merge(Class_genus, Trait_list_Class, by = "Class", all.x = TRUE)
  Class_genus$Genus <- Class_genus$Species
  Class_genus <- Class_genus[,names(Trait_genus)]
  Phylum_genus <- merge(Phylum_genus, Trait_list_Phylum, by = "Phylum", all.x = TRUE)
  Phylum_genus$Genus <- Phylum_genus$Species
  Phylum_genus <- Phylum_genus[,names(Trait_genus)]
  
  Trait_platform <- rbind(Trait_genus, Family_genus, Order_genus, Class_genus,Phylum_genus)
  
  # save amount of phyla for later use in the output file
  Phyla<-Trait_platform %>%
    group_by(Sample) %>%
    summarise(numberPhyla = n_distinct(Phylum))
  
  Phyla$TotalPhyla <- n_distinct(Trait_platform$Phylum)
  Phyla$TotalSpecies <- n_distinct(Trait_platform$Species)
  
  
  Trait_platform <- Trait_platform[,c(-1,-5:-9)]
  
  
  #######################################################################################################
  # Trait analysis
  #######################################################################################################
  
  #=============================================
  # Community weighted mean
  #=============================================
  Abundance_matrix <- Trait_platform[,c("Sample", "Species", "Abundance")]
  
  trait_matrix <- ddply(Trait_platform[,-1], .(Species), colwise(mean))
  row.names(trait_matrix) <- trait_matrix$Species
  trait_matrix <- trait_matrix[,-1:-2]
  
  Abundance_matrix_cast <- acast(Abundance_matrix, Sample ~ Species , value.var='Abundance', fun.aggregate=sum, margins=FALSE)
  Abundance_matrix_cast <- data.matrix(Abundance_matrix_cast)
  
  resCWM <- functcomp(trait_matrix, Abundance_matrix_cast, CWM.type = "all")
  
  resCWM$Sample <- row.names(resCWM)
  
  CWM_frame_average <- rowMeans(resCWM[,-ncol(resCWM)], na.rm=TRUE)
  
  CWM_cal <- as.data.frame(CWM_frame_average)
  CWM_cal$Station <- row.names(CWM_cal)
  names(CWM_cal) <- c("CWM", "Sample")
  
  #=============================================
  # making distance matrix https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
  #=============================================
  
  sr <- pcoa(gowdis(trait_matrix[, 1:6]) / max(gowdis(trait_matrix[, 1:6])))$vectors[,1]
  Body <- pcoa(gowdis(trait_matrix[, 7:12]) / max(gowdis(trait_matrix[, 7:12])))$vectors[,1]
  Size <- pcoa(gowdis(trait_matrix[, 13:16]) / max(gowdis(trait_matrix[, 13:16])))$vectors[,1]
  Reproduction <- pcoa(gowdis(trait_matrix[, 17:20]) / max(gowdis(trait_matrix[, 17:20])))$vectors[,1]
  zone <- pcoa(gowdis(trait_matrix[, 21:23]) / max(gowdis(trait_matrix[, 21:23])))$vectors[,1]
  habit <- pcoa(gowdis(trait_matrix[, 24:29]) / max(gowdis(trait_matrix[, 24:29])))$vectors[,1]
  layer <- pcoa(gowdis(trait_matrix[, 30:33]) / max(gowdis(trait_matrix[, 30:33])))$vectors[,1]
  Trophic <- pcoa(gowdis(trait_matrix[, 34:39]) / max(gowdis(trait_matrix[, 34:39])))$vectors[,1]
  Mobility <- pcoa(gowdis(trait_matrix[, 40:43]) / max(gowdis(trait_matrix[, 40:43])))$vectors[,1]
  Turbation <- pcoa(gowdis(trait_matrix[, 44:48]) / max(gowdis(trait_matrix[, 44:48])))$vectors[,1]
  
  all.dist <- cbind(data.frame(sr),data.frame(Body), data.frame(Size), data.frame(Reproduction),
                    data.frame(zone), data.frame(habit), data.frame(layer), data.frame(Trophic),
                    data.frame(Mobility), data.frame(Turbation))
  
  
  
  #=============================================
  # Testing fundiversity https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
  #=============================================
  
  #library(fundiversity)
  library(FD)
  
  resFD.alltraits <- suppressWarnings(dbFD(all.dist, Abundance_matrix_cast,
                                           calc.CWM = FALSE, stand.FRic = TRUE, messages = F, m = "min", calc.FRic = TRUE, calc.FDiv = FALSE))
  
  
  # #Scores for FRic
  FRic_frame <- data.frame(resFD.alltraits$FRic)
  FRic_frame$Sample <- row.names(FRic_frame)
  
  # Calculatin the rest of the functional indicies
  
  FDis <- fd_fdis(all.dist, Abundance_matrix_cast)
  FDiv <- fd_fdiv(all.dist, Abundance_matrix_cast)
  FEve <- fd_feve(all.dist, Abundance_matrix_cast)
  RaoQ <- fd_raoq(all.dist, Abundance_matrix_cast)
  
  
  Functional_trait <- cbind(FDis, FDiv[,2], FEve[,2], FRic_frame[,1], RaoQ[,2])
  names(Functional_trait) <- c("Sample", "FDis", "FDiv", "FEve", "FRic", "RaoQ")
  
  Functional_trait <- merge(Functional_trait, CWM_cal, by = "Sample", all = TRUE)
  
  
  print("2_BioQ - Done")
  
  
  
  
  #========================================================================
  #
  # INTERACTIONS INDICIE CALCULATIONS
  #
  # by Ole B. Brodnicke, DHI Aug 2023
  #
  #========================================================================
  
  
  #===========================================================
  # CO-OCCURENCE USING IGRAPH ACROSS STATIONS
  #===========================================================
  
  print(Species_Summary)
  station_matrix <- Species_Summary

  co_occurence <- data.frame()
  result_df<-data.frame()
  for (R in unique(Species_Summary$Sample)) {
    
    station_year_data <- station_matrix[station_matrix$Sample == R & station_matrix$Abundance !=0, ] #Removing 0's for the further analysis
    
    # Calculate dissimilarity matrix
    demo_matrix <- vegdist(station_year_data[, 3], method = "bray") #Maybe we can look at other distances
    
    # Convert dissimilarity matrix to a distance matrix
    distance_matrix <- as.dist(demo_matrix)  #The values stay the same here ? Not sure this is needed.
    
    # Create a network object
    network <- graph.adjacency(as.matrix(distance_matrix), weighted = TRUE, mode = "undirected", diag = FALSE)
    
    # Set node names
    V(network)$name <- station_year_data[, 2]
    
    # Calculate network density
    E.density <- edge_density(network, loops = T) 
    
    # Calculate diameter
    N.diam <- 1/diameter(network, directed = FALSE) #Because a smaller diameter is better we change the range so that a larger number is better (inverse)
    
    # Calculate average path length
    A.path <- 1/mean_distance(network, directed = FALSE) #Because a shorter path is better we change the range so that a larger number is better (inverse)
    
    
    
    # Append results to the dataframe
    result_df <- cbind(Sample=R, E.density, N.diam, A.path) # Adding the Sample back to the values
    
    co_occurence <- rbind(co_occurence, result_df)
    
  }
  
  # transform the avlues back to numerical variables
  co_occurence[,2:4] <- lapply(co_occurence[, 2:4], as.numeric)
  
  #===========================================================
  # Trophic-guild interaction
  #===========================================================
  
  #Trophic traits from CEFAS
  
  
  station_matrix$Genus <- unlist(lapply(strsplit(station_matrix$Species, " "), "[", 1))
  
  station_matrix_genus <- ddply(station_matrix, .( Sample, Genus), summarise, 
                                Abundance = sum(Abundance, na.rm = TRUE), Abundance = sum(Abundance, na.rm=TRUE))
  
  Platform_year_traits <- merge(station_matrix_genus, Trait_list, by = "Genus", all.x = TRUE)
  
  station_matrix_genus_2<-Platform_year_traits[,c("Genus", "Sample", "Abundance","f_Suspension","f_Surface_deposit","f_Subsurface_deposit","f_Scavenger","f_Predator", "f_Parasite")]
  
  
  # Create an entry per genus per trophic trait (f_)
  df_long <- station_matrix_genus_2 %>%
    pivot_longer(
      cols = starts_with("f_"),  # Specify the columns to pivot
      names_to = "Trophic",      # Name of the new categorical variable
      values_to = "Trait_Value"     # Name of the value column
    ) %>%
    filter(Trait_Value > 0)
  
  #Remove a column we dont need
  df_long$Trait_Value<-NULL
  
  
  combined_data <- data.frame()
  # Loop through each station and year and make unique trait combinations
  for (station in unique(df_long$Sample)) {
    station_year_data <- df_long[df_long$Sample == station, ]
    
    if (length(unique(station_year_data$Trophic)) > 1) {
      trait_combinations <- unique(t(combn(unique(station_year_data$Trophic), 2)))
      
      trait_combinations_df <- data.frame(Trait1 = trait_combinations[, 1],
                                          Trait2 = trait_combinations[, 2])
    } else { # in case there is only 1 trophic level per station:
      trait_combinations <- 0
      
      trait_combinations_df <- data.frame(Trait1 = 0,
                                          Trait2 = 0)
    }
    
    # trait_combinations <- unique(t(combn(unique(station_year_data$Trophic), 2)))}
    
    # trait_combinations_df <- data.frame(Trait1 = trait_combinations[, 1],
    #                                     Trait2 = trait_combinations[, 2])
    
    
    trait_combinations_df$ratio <- sapply(1:nrow(trait_combinations_df), function(x) {
      
      trait1_abundance <- sum(station_year_data$Abundance[station_year_data$Trophic == trait_combinations_df$Trait1[x]])
      
      trait2_abundance <- sum(station_year_data$Abundance[station_year_data$Trophic == trait_combinations_df$Trait2[x]])
      
      ratio <- trait1_abundance / trait2_abundance
      return(ratio)
    })
    
    Results <- cbind(Sample = station, trait_combinations_df)
    
    combined_data <- rbind(combined_data, Results)
  }
  
  
  
  # Calculate the the geometric mean of ratios
  agg_ratio <- combined_data %>%
    group_by(Sample) %>%
    summarize(mean_ratio = exp(mean(log(ratio))))
  
  #Merging with interactions
  Interactiv_indicies <- merge(co_occurence, agg_ratio, by = "Sample", all = TRUE)
  
  
  print("3_BioQ - Done")
  

  
  
  
  #========================================================================
  #
  # STANDARDIZATION OF INDICIE
  #
  # by Wong Xin Huei, DHI June 2023
  #
  #========================================================================
  
  #=============================================
  # Reading data from previous scripts
  #=============================================
  
  Tax_data <- Taxonomic_indicies
  Fun_data <- Functional_trait
  Inc_data <- Interactiv_indicies
  
  #Removing S and N from Tax
  # Tax_data <- subset(Tax_data, select = -c(S, NIS,n.rare , a.rare , N)) #Remove NIS n.rare etc as well?
  Tax_data <- subset(Tax_data, select = -c(S, N)) #Remove NIS n.rare etc as well?
  
  Indicie_all <- merge(Tax_data, Fun_data, by = c("Sample"))
  Indicie_all <- merge(Indicie_all, Inc_data, by = c("Sample"))
  
  Reference_site <- c("Ref.1", "Ref.2", "Ref", "Ref.") 
  
  #=============================================
  # Calculating standard site if multiple refs
  #=============================================
  
  Reference_globe <- subset(Indicie_all, Sample %in% Reference_site)
  
  Reference_globe$Sample <- "Ref_globe"
  Reference_globe <- ddply(Reference_globe, .(Sample), colwise(mean, na.rm=TRUE))
  
  #=============================================
  # Standardizing and calculating
  #=============================================
  
  Score_year <- data.frame()
  
  Score_each<-Indicie_all
  var_cal <- Indicie_all[,-1]
  var_overall <- ddply(var_cal, .(), colwise(sd, na.rm=TRUE)) # calculation of the column-wise standard deviation.
  
  #Normalizing scores
  Normalized_station <- data.frame()
  for (station in unique(Score_each$Sample)){
    Impact <- subset(Score_each, Sample == station)
    
    Results <- Impact[,"Sample"]
    for (s in names(Impact[,-1])){
      Indicie <- Impact[,s]
      Ref_val <- Reference_globe[,s]
      Score <- ifelse(var_overall[,s] > 0,(Indicie- Ref_val)/var_overall[,s], (Indicie- Ref_val)/Ref_val)
      
      Results <- cbind(Results, Score) 
      
    }
    colnames(Results) <- c("Sample", names(Impact[,-1]))
    Normalized_station <- rbind(Normalized_station, Results)
    
  }
  
  Score_year <- rbind(Score_year, Normalized_station)
  
  #Change all the numeric varibales back to numeric
  
  Score_year[,-1] <- lapply(Score_year[, -1], as.numeric)
  
  
  
  #=============================================
  # Adding BioQ calculations
  #=============================================
  
  Tax_indicies <- names(subset(Tax_data, select = -c(Sample)))
  Fun_indicies <- names(subset(Fun_data, select = -c(Sample)))
  Inc_indicies <- names(subset(Inc_data, select = -c(Sample)))
  
  BioQ_scores <- data.frame()
  
  
  Tax_Scores <-  Score_year[,c(Tax_indicies)]
  Tax_Scores[is.na(Tax_Scores)] <- 0
  Fun_Scores <-  Score_year[,c(Fun_indicies)]
  Fun_Scores[is.na(Fun_Scores)] <- 0
  Inc_Scores <-  Score_year[,c(Inc_indicies)]
  Inc_Scores[is.na(Inc_Scores)] <- 0
  
  Tax_Score <- rowMeans(Tax_Scores)
  Fun_Score <- rowMeans(Fun_Scores)
  Inc_Score <- rowMeans(Inc_Scores)
  
  Score_year$Tax_mean <- Tax_Score
  Score_year$Fun_mean <- Fun_Score
  Score_year$Inc_mean <- Inc_Score
  
  
  Score_year$BioQ <-  rowMeans(cbind(Tax_Score, Fun_Score, Inc_Score))
  
  BioQ_scores<-Score_year
  
  
  # Gettign S and N for the visualisation:
  
  Abundance <- Taxonomic_indicies
  Abundance <- subset(Abundance, select = c(Sample, S, N))
  
  Abundance <- merge(Abundance, Phyla, by = c("Sample"))
  
  Abundance$TotalN <- sum(Abundance$N)
  
  BioQ_scores <- merge(BioQ_scores, Abundance, by = c("Sample"))
  
  #write.csv(BioQ_scores, "./Indput-Output/Standardized_BioQ.csv")
  
  print("4_BioQ - Done")
  # print("BioQ_scores calculated")
  
  
  
  # Perform operations on the data
  summary_data <- BioQ_scores
  
  # Return the summary data
  return(BioQ_scores)
  
  
}


