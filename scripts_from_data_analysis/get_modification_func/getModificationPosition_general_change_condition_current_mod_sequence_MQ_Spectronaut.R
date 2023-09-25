
# mod_seq <- "_S[Phospho (STY)]GGQRHSPLSQR_.3"
# mod_seq <- "_S[Phospho (STY)]SSFREM[Oxidation (M)]DGQPER_.3"

getModificationPosition_general <- function( mod_seq, software_name){ 
  if (software_name == "Spectronaut"){
    split_mod_seq <- strsplit(mod_seq, ".", fixed = TRUE)[[1]]
    
    
    mod_seq1 <- substr(split_mod_seq[1], 2, nchar(split_mod_seq[1])-1)
    mod_seq1 <- gsub("\\[", "(", mod_seq1)
    mod_seq1 <- gsub("\\]", ")", mod_seq1)
    
    modification_labels <- base::regmatches(mod_seq1, gregexpr("\\(.*?\\)\\)", mod_seq1))[[1]]
    naked_peptide <- gsub("\\(.*?\\)|\\)", "", mod_seq1)
    
    
    # modification_labels <- base::regmatches(mod_seq1, gregexpr("\\[.*?\\]", mod_seq1))[[1]]
    # naked_peptide <- gsub("\\[.*?\\]", "", mod_seq1)
    # 
    # 
    
    modification_index_list <- list()
    modification_index_list$pep_seq <- naked_peptide
    modification_index_list$naked_peptide_length <- nchar(naked_peptide)
    modification_index_list$charge <- split_mod_seq[2]
    
  } else if (software_name == "MQ_214"){
    
    
    mod_seq1 <- substr(mod_seq, 2, nchar(mod_seq)-1)
    
    modification_labels <- base::regmatches(mod_seq1, gregexpr("\\(.*?\\)\\)", mod_seq1))[[1]]
    naked_peptide <- gsub("\\(.*?\\)|\\)", "", mod_seq1)
    
    modification_index_list <- list()
    modification_index_list$pep_seq <- naked_peptide
    modification_index_list$naked_peptide_length <- nchar(naked_peptide)
    
  }
  else{
    print("Your selection could not find in this function. Make sure that either selecting MQ_214 and Spectronaut")
  }
  
  for ( mod in unique(modification_labels) ) {
    
    ## ## Remove other modifications from the sequence 
    ## that are not being accessed to get a more accurate position index
    
    mods_to_remove_from_sequence <- modification_labels[ !(modification_labels %in% mod) ]
    
    ## In case of 2 PTMs allowed
    if( length(mods_to_remove_from_sequence)==1 ){
      
      current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence )), "", mod_seq1 )
    }
    ## In case of 3 PTMs allowed (should be max value by default)
    else if (length(mods_to_remove_from_sequence) == 2) {
      mid_mod_sequence <- gsub(gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence[1])), "", mod_seq1)
      current_mod_sequence <- gsub(gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence[2])), "", mid_mod_sequence)
    }
    else {
      #current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence )), "", mod_seq1 )
      current_mod_sequence <-mod_seq1
    } 
    ## Replace current Modification string with an identifier to get amino acid index
    current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mod )), "@", current_mod_sequence )
    
    ## Get positions of current modification based on identifier
    pos_mod <- gregexpr( '@', current_mod_sequence )
    
    ## Get actual amino acid index
    if ( as.numeric(unlist((pos_mod)))[1] == 1 ){
      pos_mod <- pos_mod
      
    }else if (length(as.numeric(unlist((pos_mod))))==1){
      pos_mod <-as.numeric(unlist((pos_mod))) -1
      
    }else if (length(as.numeric(unlist((pos_mod))))==2){
      
      pos_mod <-as.numeric(unlist((pos_mod))) -1
      pos_mod[2] <- pos_mod[2] -1
      
      
    }else if(length(as.numeric(unlist((pos_mod))))==3){
      
      pos_mod <- as.numeric(unlist((pos_mod))) -1
      
      pos_mod[2] <- pos_mod[2] -1
      pos_mod[3] <- pos_mod[3] -2
    }
    
     
    
    ## Store position in list
    modification_index_list[[paste0('modification_',mod)]] <- as.numeric(unlist(pos_mod))
    
  }
  return( modification_index_list )    
  
}