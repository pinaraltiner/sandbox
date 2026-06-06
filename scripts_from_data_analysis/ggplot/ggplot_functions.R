#
library(ggplot2)
#library(dplyr)
library(gghalves)
###### SPECTRONAUT & POSITION EXTRACTION ######
################
find_max_value_and_pos <- function(ptm_count, ptm_prob, ptm_pos) {
  
  prob_values <- as.numeric(unlist(strsplit(ptm_prob, ";")))
  pos_values <- as.numeric(unlist(strsplit(ptm_pos, ";")))
  result <- list()
  dim_count <- dim(as.data.frame(str_match_all(pattern = "\\[Phospho", ptm_count)))[1]
  max_prob <- max(prob_values)
  max_pos <- pos_values[which(prob_values == max_prob)]
  
  
  if (dim_count >= 2) {
    
    index_order <- order(prob_values, decreasing = TRUE)
    
    # Sort prob_values and pos_values using the same index order
    sorted_prob_values <- prob_values[index_order]
    sorted_pos_values <- pos_values[index_order]
    
    highest_dim <- sorted_prob_values[1:dim_count]
    
    #for (i in 1:length(highest_dim)){
    result["prob"] <- max(sorted_prob_values)#paste(sorted_prob_values[1:length(highest_dim)],collapse ="&")
    result["pos"] <- paste(sorted_pos_values[1:length(highest_dim)],collapse ="&")
    result$mod <- "two_phospho"
    #result<-list.append(paste(max_prob, max_pos, sep = "_"))
    #result$position <- paste(max_prob, max_pos, sep = "_")
    #result$position1 <- paste(max_prob, max_pos, sep = "_")
    #positions = paste0("position",i)
    #probs=paste0("prob",i)
    #max_pos <- 
    #result[positions] <- paste(pos_values[which(prob_values[i+1] == highest_dim[i])],pos_values[which(prob_values[i+1] == highest_dim[i])],sep = "&")
    #append(result[[positions]], paste(pos_values[i],pos_values[i],sep = "&"))
    #result[[probs]] <- append(result[[probs]], paste(sorted_probs[i],sorted_probs[i],sep = "&"))
    #}
  }else if (length(max_pos) == 1) {
    result$prob <- max_prob
    result$pos <- max_pos
    result$mod <- "mono phospho"
  } else {
    #### IF THIS PART CREATES AN ERROR, REMOVE THE POS AND 
    #### SCORE VALUE JUST RETURNED " Non-Distinguishable"
    for (i in 1:length(max_pos)){
      #result<-list.append(paste(max_prob, max_pos, sep = "_"))
      #result$position <- paste(max_prob, max_pos, sep = "_")
      #result$position1 <- paste(max_prob, max_pos, sep = "_")
      probs = paste0("prob",i)
      poses = paste0("pos",i)
      result[[probs]] <- append(result[[probs]], paste(max_pos[i]))
      result[[poses]] <- append(result[[poses]], paste(max_prob[i]))
      
      result$mod <- "non-distinguishable"
    }
  } 
  return(result)
}
################


###### Proteome Discoverer MODIFICATION & POSITION EXTRACTION ######
################
extract_phospho_numbers <- function(input_string) {
  phospho_part <- regmatches(input_string, gregexpr("Phospho \\[[^]]+\\]", input_string))
  
  if (length(phospho_part) == 0) {
    return(list(NULL, NULL))
  }
  
  phospho_text <- phospho_part[[1]]
  letter_numbers <- gregexpr("[A-Z](\\d+)", phospho_text)
  extracted_letters <- regmatches(phospho_text, letter_numbers)[[1]]
  extracted_values <- gregexpr("\\((\\d+(?:\\.\\d+)?)\\)", phospho_text)
  extracted_values <- regmatches(phospho_text, extracted_values)[[1]]
  
  extracted_values <- gsub("\\(|\\)", "", extracted_values)  # Remove parentheses
  extracted_letters <- gsub("[A-Z]", "", extracted_letters)  # Remove letters
  ## Selecting the max is necessary for applying 
  ##           filtering at localization score
  combined_results <- list(paste(extracted_letters, collapse = "&"), paste(max(extracted_values), collapse = "&"))
  return(combined_results)
}

###### PROLINE MODIFICATION & POSITION EXTRACTION ###### 

proline_phospho_pos_extraction <- function(df){
  # Extraction of phospho positions from quant peptides object
  phospho_ptm_pos <- lapply(df, function(each_ptm_protein_positions) {
    
    ptm_list <- as.list(strsplit(each_ptm_protein_positions,"; ", fixed=TRUE)[[1]]) # Split ptm_protein_position depending on ";"
    ptm_list <- ptm_list[grepl("Phospho", ptm_list, fixed = TRUE)] # Extract only which contains "Phospho"
    phospho_positions <- lapply(ptm_list, function(ptm) { 
      #sub('Phospho \\(([A-Z]\\d+)\\)', "\\d+", ptm) #then, remove "Phospho" and remain only positions
      as.character(str_extract(ptm, "\\d+"))
      
    })
    
    phospho_positions_as_str <- paste(phospho_positions, collapse="&") #combine each position with "|"
    
  })
  # Data conversion 
  phospho_ptm_pos_df <- t(as.data.frame(phospho_ptm_pos))
  rownames(phospho_ptm_pos_df) <- 1:length(phospho_ptm_pos_df)
  
  return(phospho_ptm_pos_df)
  
}



##Distribution of both mean abundance 
#                                    and Exp. Ratios 
#                                                   of every sample with two pools

library(dplyr)
### DATA FILTERING ###
filter_NA <- function(df, samp_names, num_allowed_NA, num_expected_nonNA) {
  # Create a logical matrix indicating non-NA values for the specified columns
  nonNA_matrix <- sapply(samp_names, function(samp_name) {
    rowSums(!is.na(select(df, contains(samp_name)))) >= num_expected_nonNA
  })
  
  # Convert the logical matrix to a data frame
  nonNA_df <- as.data.frame(nonNA_matrix)
  
  # Filter rows where the number of TRUE values is >= num_allowed_NA
  valid_rows <- rowSums(nonNA_df) >= num_allowed_NA
  
  # Return the filtered data frame
  filtered_df <- df[valid_rows, ]
  
  return(filtered_df)
}




### CONDITIONAL IMPUTATION ###
conditional_imputation <- function(df,impute_val,samp_names,num_NA_imputed){
  all_complete_df <- NULL
  for (i in 1:length(samp_names)){
    assign(paste0("df_",samp_names[i]), df %>%
             mutate(impute=ifelse(rowSums(is.na(select(df,contains(samp_names[i]))))>=num_NA_imputed,TRUE,FALSE)) %>%
             select(contains(samp_names[i]) | contains("impute"),indx))
    
    sel_for_impt <- get(paste0("df_",samp_names[i])) %>%
      filter(grepl(TRUE,impute)) %>%
      select(!impute)
    
    nonimp_df <- get(paste0("df_",samp_names[i])) %>%
      filter(!grepl(TRUE,impute)) %>%
      select(!impute)
    
    sel_impute_values <- impute_val[(3*(i-1)+1):(3*(i-1)+3)]
    
    for(k in 1:length(sel_impute_values)){
      sel_for_impt[,k]<-sel_impute_values[k]
    }
    
    assign(paste0("df_complete",samp_names[i]),bind_rows(nonimp_df, sel_for_impt))
    assign(paste0("ordered_df_complete",i),get(paste0("df_complete",samp_names[i]))[order(get(paste0("df_complete",samp_names[i]))$indx),])
    
    all_complete_df <- all_complete_df %>% bind_cols(get(paste0("ordered_df_complete",i))) %>% select(!indx)
    rm(sel_for_impt,nonimp_df)
  }
  return(all_complete_df)
}
############## ############## ############## ##############
## FUNCTION FOR CV CALCULATION ##### FIRST VERSION ##############
############## ############## ############## ##############
# CV_calculator <- function(abundance_col,
#                           data,
#                           replicate,
#                           iter,
#                           exp_id,
#                           acquisiton_type,software_name
#                           
# ){
#   for (i in 1:iter){ ## it has to be started from 2 for experiment 3
#     
#     old_mean_name <- paste0(abundance_col,i)
#     new_name <- paste0("CV",i)
#     
#     if(exp_id ==3){
#       dftmp_mean <- data %>%
#         select(pep_with_pos,species,contains(paste0(abundance_col,i))) %>% 
#         rename_with(~"row_mean", contains(abundance_col)) %>%
#         #bind_cols(1:nrow(data)) %>%
#         drop_na() 
#       #rename(index = (replicate+1)) 
#       
#       if (i==2){
#         
#         assign(paste0('A_CV',i), data %>%
#                  select(contains(paste0("E3-A",i))) %>%
#                  drop_na() %>%
#                  #rowwise() %>%
#                  mutate(
#                    row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
#                    # Add other transformations here
#                  ) %>%
#                  #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
#                  bind_cols(dftmp_mean) %>%
#                  mutate(CV=(row_sd/row_mean)*100) %>%
#                  bind_cols(1:nrow(dftmp_mean)) %>%
#                  rename(index = ((replicate)+(2*(replicate)))) %>%
#                  #rename_with(~paste0(new_name,.x,recycle0 = FALSE),starts_with("CV")) %>%
#                  rename_with(~new_name, contains("CV")) %>%
#                  select(!contains("E3-A")))
#       }else{
#         
#         assign(paste0('A_CV',i), data %>%
#                  select(contains(paste0("E3-A",i))) %>%
#                  drop_na() %>%
#                  #rowwise() %>%
#                  mutate(
#                    row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
#                    # Add other transformations here
#                  ) %>%
#                  #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
#                  bind_cols(dftmp_mean) %>%
#                  mutate(CV=(row_sd/row_mean)*100) %>%
#                  bind_cols(1:nrow(dftmp_mean)) %>%
#                  rename(index = (replicate+6)) %>%
#                  select(CV,index) %>%
#                  rename_with(~new_name, contains("CV")))
#       }
#       print(i)
#       
#     }else if (exp_id==2){
#       dftmp_mean <- data %>%
#         select(pep_with_pos,Pool,contains(paste0(abundance_col,i))) %>% 
#         rename_with(~"row_mean", contains(abundance_col)) %>%
#         #bind_cols(1:nrow(data)) %>%
#         drop_na() 
#       #rename(index = (replicate+1)) 
#       
#       if (i==1){
#         
#         assign(paste0('A_CV',i), data %>%
#                  select(contains(paste0("E2-A",i))) %>%
#                  drop_na() %>%
#                  #rowwise() %>%
#                  mutate(
#                    row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
#                    # Add other transformations here
#                  ) %>%
#                  #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
#                  bind_cols(dftmp_mean) %>%
#                  mutate(CV=(row_sd/row_mean)*100) %>%
#                  bind_cols(1:nrow(dftmp_mean)) %>%
#                  rename(index = (replicate+(2*(replicate)))) %>%
#                  #rename_with(~paste0(new_name,.x,recycle0 = FALSE),starts_with("CV")) %>%
#                  rename_with(~new_name, contains("CV")) %>%
#                  select(!contains("E2-A")))
#       }else{
#         
#         assign(paste0('A_CV',i), data %>%
#                  select(contains(paste0("E2-A",i))) %>%
#                  drop_na() %>%
#                  #rowwise() %>%
#                  mutate(
#                    row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
#                    # Add other transformations here
#                  ) %>%
#                  #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
#                  bind_cols(dftmp_mean) %>%
#                  mutate(CV=(row_sd/row_mean)*100) %>%
#                  bind_cols(1:nrow(dftmp_mean)) %>%
#                  rename(index = (replicate+6)) %>%
#                  select(CV,index) %>%
#                  rename_with(~new_name, contains("CV")))
#       }
#       print(i)
#     }
#   
#   }
#   
#   df_CV <- A_CV2 %>% #A_CV1 %>%
#     select(!row_sd) %>%
#     #left_join(A_CV2,by="index") %>%
#     left_join(A_CV3,by="index") %>%
#     left_join(A_CV4,by="index") %>%
#     left_join(A_CV5,by="index") %>%
#     pivot_longer(cols = contains("CV"),
#                  names_to = "CV_samples",
#                  values_to = "CV_values") %>%
#     #filter(!grepl("Unexpected",Pool)) %>%
#     #mutate(CV_values=(CV_values*100)) %>%
#     mutate(acq_type=acquisiton_type) %>%
#     mutate(soft_name=software_name) %>%
#     drop_na()
#   
#   return(df_CV)
# }

############## ############## ############## ##############
## FUNCTION FOR CV CALCULATION ##### SECOND VERSION ##############
############## ############## ############## ##############
CV_calculator <- function(abundance_col,
                          data,
                          replicate,
                          iter,
                          exp_id,
                          acquisiton_type,software_name
                          
){
  
  
  
  if(exp_id ==3){
    
    for(i in 2:iter){
      
      old_mean_name <- paste0(abundance_col,i)
      new_name <- paste0("CV",i)
      
      dftmp_mean <- data %>%
        select(pep_with_pos,species,contains(paste0(abundance_col,i))) %>% 
        rename_with(~"row_mean", contains(abundance_col)) %>%
        #bind_cols(1:nrow(data)) %>%
        drop_na() 
      #rename(index = (replicate+1)) 
      
      if (i==2){
        
        assign(paste0('A_CV',i), data %>%
                 select(contains(paste0("E3-A",i))) %>%
                 drop_na() %>%
                 #rowwise() %>%
                 mutate(
                   row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
                   # Add other transformations here
                 ) %>%
                 #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
                 bind_cols(dftmp_mean) %>%
                 mutate(CV=(row_sd/row_mean)*100) %>%
                 bind_cols(1:nrow(dftmp_mean)) %>%
                 rename(index = ((replicate)+(2*(replicate)))) %>%
                 #rename_with(~paste0(new_name,.x,recycle0 = FALSE),starts_with("CV")) %>%
                 rename_with(~new_name, contains("CV")) %>%
                 select(!contains("E3-A")))
      }else{
        
        assign(paste0('A_CV',i), data %>%
                 select(contains(paste0("E3-A",i))) %>%
                 drop_na() %>%
                 #rowwise() %>%
                 mutate(
                   row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
                   # Add other transformations here
                 ) %>%
                 #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
                 bind_cols(dftmp_mean) %>%
                 mutate(CV=(row_sd/row_mean)*100) %>%
                 bind_cols(1:nrow(dftmp_mean)) %>%
                 rename(index = (replicate+6)) %>%
                 select(CV,index) %>%
                 rename_with(~new_name, contains("CV")))
      }
      print(i)
    }
    
    df_CV <- A_CV2 %>% #A_CV1 %>%
      select(!row_sd) %>%
      #left_join(A_CV2,by="index") %>%
      left_join(A_CV3,by="index") %>%
      left_join(A_CV4,by="index") %>%
      left_join(A_CV5,by="index") %>%
      pivot_longer(cols = contains("CV"),
                   names_to = "CV_samples",
                   values_to = "CV_values") %>%
      #filter(!grepl("Unexpected",Pool)) %>%
      #mutate(CV_values=(CV_values*100)) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>%
      drop_na()
    
    
  }else if (exp_id==2){
    
    for (i in 1:iter){
      
      old_mean_name <- paste0(abundance_col,i)
      new_name <- paste0("CV",i)
      
      dftmp_mean <- data %>%
        select(pep_with_pos,Pool,contains(paste0(abundance_col,i))) %>% 
        rename_with(~"row_mean", contains(abundance_col)) %>%
        #bind_cols(1:nrow(data)) %>%
        drop_na() 
      #rename(index = (replicate+1)) 
      
      if (i==1){
        
        assign(paste0('A_CV',i), data %>%
                 select(contains(paste0("E2-A",i))) %>%
                 drop_na() %>%
                 #rowwise() %>%
                 mutate(
                   row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
                   # Add other transformations here
                 ) %>%
                 #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
                 bind_cols(dftmp_mean) %>%
                 mutate(CV=(row_sd/row_mean)*100) %>%
                 bind_cols(1:nrow(dftmp_mean)) %>%
                 rename(index = (replicate+(2*(replicate)))) %>%
                 #rename_with(~paste0(new_name,.x,recycle0 = FALSE),starts_with("CV")) %>%
                 rename_with(~new_name, contains("CV")) %>%
                 select(!contains("E2-A")))
      }else{
        
        assign(paste0('A_CV',i), data %>%
                 select(contains(paste0("E2-A",i))) %>%
                 drop_na() %>%
                 #rowwise() %>%
                 mutate(
                   row_sd = apply(select(., where(is.numeric)), 1, sd, na.rm = TRUE),
                   # Add other transformations here
                 ) %>%
                 #mutate(row_sd = sd(c_across(where(is.numeric)), na.rm = TRUE)) %>%
                 bind_cols(dftmp_mean) %>%
                 mutate(CV=(row_sd/row_mean)*100) %>%
                 bind_cols(1:nrow(dftmp_mean)) %>%
                 rename(index = (replicate+6)) %>%
                 select(CV,index) %>%
                 rename_with(~new_name, contains("CV")))
      }
      print(i)
      
    }
    
    df_CV <- A_CV1 %>% 
      select(!row_sd) %>%
      left_join(A_CV2,by="index") %>%
      left_join(A_CV3,by="index") %>%
      left_join(A_CV4,by="index") %>%
      left_join(A_CV5,by="index") %>%
      pivot_longer(cols = contains("CV"),
                   names_to = "CV_samples",
                   values_to = "CV_values") %>%
      #filter(!grepl("Unexpected",Pool)) %>%
      #mutate(CV_values=(CV_values*100)) %>%
      mutate(acq_type=acquisiton_type) %>%
      mutate(soft_name=software_name) %>%
      drop_na()
    
  }
  
  
  return(df_CV)
}


##################################
### ADDING NEUTRAL MASS FOR PEPTIDE TO MO MAPPING ###

# calculate_sum <- function(sequence, positions, mapping_df) {
#   # Split the sequence into individual letters
#   letters <- unlist(strsplit(sequence, ""))
#   
#   # Map the letters to their corresponding numbers
#   values <- mapping_df$mono_isotopic[match(letters, mapping_df$letter)]
#   
#   # Calculate the base sum of the values
#   base_sum <- sum(values, na.rm = TRUE)
#   
#   # Check if the positions column contains "&" or only numeric values
#   if (grepl("&", positions)) {
#     additional_value <- 177.9568
#   } else {
#     additional_value <- 97.97686
#   }
#   
#   # Add the additional value based on the condition
#   final_sum <- base_sum + additional_value
#   
#   return(final_sum)
# }

calculate_sum <- function(Sequence, Positions, mapping_df) {
  # Split the sequence into individual letters
  letters <- unlist(strsplit(Sequence, ""))
  
  # Map the letters to their corresponding numbers
  values <- mapping_df$mono_isotopic[match(letters, mapping_df$letter)]
  
  # Calculate the base sum of the values
  base_sum <- sum(values, na.rm = TRUE) +18
  
  # Check if the positions column contains "&" or only numeric values
  if (grepl("&", Positions)) {
       # Count the number of "&" symbols in the positions string
       count_ampersand <- nchar(gsub("[^&]", "", Positions))
       additional_value <-(count_ampersand+1) * 80
  }else {
       additional_value <- 80
       }
  # Add the additional value to the base sum
  final_sum <- base_sum + additional_value
  
  return(final_sum)
}

##################################

###################
gg_quant_ratio_acc <- function(data_set,
                               color_df,
                               x_df,
                               y_df,
                               median_col,
                               xmin,
                               header,
                               x_lab,
                               y_lab,
                               color_lab,
                               subtitle_txt
                               ){
  
  ggplot(data_set,aes(x=x_df,y=y_df,color=color_df,xmin)) + 
    geom_point(size=12) + stat_summary(fun.y=median, geom="point", shape=18,
                              size=15, color=median_col) + 
    geom_smooth(method = "lm",col = "black",size=1.9) +
    theme_minimal() +
    
    theme(legend.text = element_text(size = 45),
          axis.title.x = element_text(size = 45),
          axis.title.y = element_text(size = 45),
          plot.title = element_text(size = 55),
          legend.title = element_text(size = 45),
          axis.text.x = element_text(size = 45),
          axis.title = element_text(size = 45),
          axis.text.y = element_text(size = 45),
          plot.subtitle = element_text(size = 45),
          strip.text.x = element_text(
            size = 20))+
    ggtitle(header) +
    labs(x=x_lab, y=y_lab,color= color_lab,subtitle = subtitle_txt) +
    scale_color_brewer(palette = 1,direction=-1) #+
    # scale_y_continuous(
    #   limits = c(-5,9), 
    #   breaks = seq(-5, 9,1)
    # ) +
    # scale_x_continuous(
    #   limits = c(xmin,7), 
    #   breaks = seq(xmin, 7,0.5)
    # ) 
}

gg_density <- function(data_set,
                                 x_df,
                                 fill_df,
                                 color_df,
                                 header,
                                 facet_df,
                                 x_lab,
                                 fill_lab,
                                 color_lab,
                       subtitle_txt){
  data_set$facet <- data_set[[facet_df]]
  ggplot(data_set,aes(x=log10(x_df), fill=fill_df )) +
    stat_density(position = "stack",alpha=0.7, aes(color=color_df),size=1)+
    #stat_density(position = "identity") 
    #https://ggplot2.tidyverse.org/reference/geom_density.html
    scale_fill_brewer(palette = 1,direction=-1) +
    scale_linetype_manual(values=c("dashed", "dotted")) +
    theme_minimal() +
    theme(legend.text = element_text(size=30), 
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30),
          strip.text.x = element_text(
            size = 15
          )
          ) +
    ggtitle(header) +
    facet_wrap(~facet, nrow = 2,scales = "free_x") +
    labs(x=x_lab,fill = fill_lab, color= color_lab,subtitle = subtitle_txt)
    
  
} 

            
    


#### THIS IS NOT TESTED ####

# gg_density_mean_abun <- function(data_set,
#                                  x_df,
#                                  y_df,
#                                  color_df,
#                                  header,
#                                  facet_df,
#                                  header,
#subtitle_txt){
#   data_set$facet <- data_set[[facet_df]]
#   ggplot(data_set, aes(x =log10(x_df) , y = log10(y_df), color=color_df)) +
#   geom_point()+
#   #facet_wrap(vars(facet_df))  + 
#   geom_smooth(formula = y ~ x,method = "loess", colour = "green", fill = "green") +
#   theme_minimal() +
#   theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
#         axis.title.x = element_text(size=30),axis.title.y = element_text(size=30),
#         plot.title = element_text(size=35),
#         legend.title=element_text(size=15),
#         axis.text=element_text(size=15),
#         axis.title=element_text(size=15)
#   ) +labs(subtitle = subtitle_txt)

### BOX-PLOT: Experimental Quantity Ratio of Synthetic Peptides  

gg_boxplt_exp_ratio <- function(data_set,
                                 x_df,
                                 y_df,
                                 fill_df,
                                 header,
                                 x_lab,
                                 y_lab,
                                 fill_lab,
                                subtitle_txt){
  ggplot(data_set,aes(x =x_df , y =log2(y_df),fill = fill_df)) +
    geom_boxplot() +
    scale_y_continuous(breaks = seq(from =round(min(log2(y_df))), to=(round(max(log2(y_df)))+2),by=1)) +
    theme_minimal() +
    theme(legend.text = element_text(size=30), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)
    ) +   stat_boxplot(geom = "errorbar") + 
    ggtitle(header) +
    labs(x=x_lab,y=y_lab,fill = fill_lab,subtitle = subtitle_txt) +
    scale_fill_brewer(palette="Dark2")
  
}

library(gghalves)

### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  
 
gg_half_boxplt_exp_ratio <- function(data_set,
                                x_df,
                                y_df,
                                fill_df,
                                header,
                                x_lab,
                                y_lab,
                                fill_lab,
                                subtitle_txt
                                ){
  ggplot(data_set,aes(x =x_df , y =log2(y_df),fill = fill_df)) +
    geom_half_boxplot(outlier.shape = NA) +
    geom_half_point(alpha = 1, show.legend = FALSE, aes(color=fill_df))+
    scale_y_continuous(breaks = seq(from =round(min(log2(y_df))), to=(round(max(log2(y_df)))+2),by=1)) +
    theme_minimal() +
    theme(legend.text = element_text(size=15),
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)
    ) + scale_fill_brewer(palette = "Dark2") +
        scale_color_brewer(palette = "Dark2")+ 
    labs(x=x_lab,y=y_lab,fill=fill_lab,subtitle = subtitle_txt) +
    ggtitle(header)
}



### HALF-BOX-PLOT & HALF-SCATTER-PLOT: Experimental Quantity Ratio of Synthetic Peptides  

gg_half_boxplt_exp_ratio_nolog <- function(data_set,
                                     x_df,
                                     y_df,
                                     fill_df,
                                     header,
                                     x_lab,
                                     y_lab,
                                     fill_lab,
                                     subtitle_txt
){
  ggplot(data_set,aes(x =x_df , y =y_df,fill = fill_df)) +
    geom_half_boxplot(outlier.shape = NA) +
    geom_half_point(alpha = 1, show.legend = FALSE, aes(color=fill_df))+
    scale_y_continuous(breaks = seq(from =round(min(y_df)), to=(round(max(y_df))),by=10)) +
    theme_minimal() +
    theme(legend.text = element_text(size=15),
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)
    ) + scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2")+ 
    labs(x=x_lab,y=y_lab,fill=fill_lab,subtitle = subtitle_txt) +
    ggtitle(header)
}






### VIOLIN-PLOT: Experimental Quantity Ratio of Synthetic Peptides   

gg_violin_exp_ratio <- function(data_set,
                                     x_df,
                                     y_df,
                                     fill_df,
                                     header,
                                     x_lab,
                                     y_lab,
                                     fill_lab,
                                     trim,
                                subtitle_txt){
  
  ggplot(data_set,aes(x =x_df , y =log2(y_df),fill = fill_df)) +
    geom_violin(trim = trim) +
    scale_y_continuous(breaks = seq(from =round(min(log2(y_df))), to=(round(max(log2(y_df)))+2),by=1)) +
    theme_minimal() +
    theme(legend.text = element_text(size=30),
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)
    ) + scale_fill_brewer(palette="Dark2") +
    labs(x=x_lab,y=y_lab, fill= fill_lab,subtitle = subtitle_txt) +
    ggtitle(header)
}

gg_violin_exp_ratio_nolog <- function(data_set,
                                x_df,
                                y_df,
                                fill_df,
                                header,
                                x_lab,
                                y_lab,
                                fill_lab,
                                trim,
                                subtitle_txt){
  
  ggplot(data_set,aes(x =x_df , y =y_df,fill = fill_df)) +
    geom_violin(trim = trim) +
    scale_y_continuous(breaks = seq(from =round(min(y_df)), to=(round(max(y_df))),by=5)) +
    theme_minimal() +
    theme(legend.text = element_text(size=30),
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)
    ) + scale_fill_brewer(palette="Dark2") +
    labs(x=x_lab,y=y_lab, fill= fill_lab,subtitle = subtitle_txt) +
    ggtitle(header)
}


#################################

gg_barplt_id_pep_count <- function(data_set,
                                x_df,
                                fill_df,
                                header,
                                ymax,
                                x_lab,
                                y_lab,
                                fill_lab,
                                caption_lab, # "NA values are removed.",
                                subtitle_txt,
                                size_num){


ggplot(data_set, aes(x=x_df, fill=fill_df)) +  geom_bar(position = position_dodge2(preserve = "single")) +
  scale_fill_brewer(palette = 1,direction=-1) +
  theme_minimal() +
  theme(legend.text = element_text(size=30), 
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        plot.title = element_text(size=35),
        plot.subtitle = element_text(size = 25),
        plot.caption = element_text(size = 25),
        legend.title=element_text(size=30),
        axis.text=element_text(size=30),
        axis.title=element_text(size=30)) +
  ggtitle(header) + 
  scale_y_continuous(breaks = seq(from=0, to=ymax,by=1000)) + #,expand = expansion(mult = c(0, 0.05))
  geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.5,size=size_num) +
  labs(x=x_lab,y=y_lab, fill= fill_lab, caption = caption_lab,subtitle = subtitle_txt)
}
  
gg_barplt_id_pep_count_stack <- function(data_set,
                                   x_df,
                                   fill_df,
                                   header,
                                   ymax,
                                   x_lab,
                                   y_lab,
                                   fill_lab,
                                   caption_lab, # "NA values are removed.",
                                   subtitle_txt,
                                   size_num){
  
  
  ggplot(data_set, aes(x=x_df, fill=fill_df)) + geom_bar(position = "stack") +
    scale_fill_brewer(palette = 1,direction=-1) +
    theme_minimal() +
    theme(legend.text = element_text(size=30), 
          axis.title.x = element_text(size=30),
          axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 25),
          plot.caption = element_text(size = 25),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)) +
    ggtitle(header) + 
    scale_y_continuous(breaks = seq(from=0, to=ymax,by=500)) +
    geom_text(color="black",aes(label=after_stat(count)),stat = "count",size=size_num,position = position_stack(vjust = 0.5))+
    labs(x=x_lab,y=y_lab, fill= fill_lab, caption = caption_lab,subtitle = subtitle_txt)
}


library(ggdist)

gg_raincloud <- function(data_set,
                         x_df,
                         y_df,
                         fill_df,
                         
                         header,
                         x_lab,
                         y_lab,
                         fill_lab,
                         caption_lab,
                         subtitle_txt){
  ggplot(data_set, aes(x= factor(x_df),
                                         y = log2(y_df),
                                         fill=factor(fill_df)
                       )) +
    ggdist::stat_halfeye(side="right", adjust = 0.5,
                          justification = -0.15,
                          .width = 0,
                          point_colour=NA,
                          alpha=0.75) +
    geom_boxplot(
                 width = .24, 
                 outlier.colour = NA,
                 alpha=1) +
    # ggdist::stat_dots(side="left",
    #                   justification= 1.1,
    #                   binwidth =.025) + 
    theme_light() + 
    theme(legend.text = element_text(size=30), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size=30),axis.title.y = element_text(size=30),
          plot.title = element_text(size=35),
          plot.subtitle = element_text(size = 25),
          plot.caption = element_text(size = 25),
          legend.title=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title=element_text(size=30)) + 
    ggtitle(header) + 
    labs(x=x_lab,y=y_lab, fill= fill_lab, caption = caption_lab,subtitle = subtitle_txt) +
    #+ #,"cyan3"
    scale_fill_brewer(palette="Dark2")

  
}
  

gg_volcano <- function(data_set,
                       x_df,
                       y_df,
                       facet_df,
                       color_df,
                       header,
                       x_lab,
                       y_lab,
                       color_lab,
                       caption_lab,
                       subtitle_txt){
  data_set$facet <- data_set[[facet_df]]
  ggplot(data_set ,aes(x =x_df, y = y_df, color=color_df)) +
    geom_point(size = 1) + #, aes(shape=merge_stat_df_final$species)
    facet_wrap(~facet) +
    scale_x_continuous(limits = c(-10, 20),breaks = seq(from = -10, to = 20, by = 2)) +  # Set the ticks for the x-axis
    scale_y_continuous(limits = c(0, 10),breaks = seq(from = 0, to = 10, by = 2))+
    #scale_y_continuous(limits = c(round(min(y_df),2), round(max(y_df),2)), breaks = seq(round(min(y_df),2), round(max(y_df),2), by = 1)) +
    #scale_x_continuous(limits = c(min(x_df), max(x_df)),breaks = seq(min(x_df), max(x_df), by = 1)) +
    scale_color_brewer(palette = "Dark2") +
    #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
    theme_bw() +
    theme(legend.text = element_text(size = 30),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          plot.subtitle = element_text(size=25),
          plot.title = element_text(size = 35),
          legend.title = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.title = element_text(size = 30),
          axis.text.y = element_text(size = 30)) +
    #expand_limits(x = 0, y = 0) +
    #geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
    #geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    labs(x=x_lab, y= y_lab, title = header, color = color_lab, subtitle = subtitle_txt) +
    guides(color = guide_legend(override.aes = list(size = 10))) 
}

  
  
  
  
  