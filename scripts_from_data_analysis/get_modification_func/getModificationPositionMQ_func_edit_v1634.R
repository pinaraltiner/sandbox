getModificationPosition_MQ <- function( mod_seq){ 
  
  mod_seq_tmp <- gsub("(?<!\\(ox\\))p", "(p)", mod_seq, perl = TRUE)
  mod_seq1 <- substr(mod_seq_tmp, 2, nchar(mod_seq_tmp)-1)
  
  modification_labels <- base::regmatches(mod_seq1, gregexpr("\\([^()]+\\)", mod_seq1))[[1]]
  naked_peptide <- gsub("\\(.*?\\)|\\)", "", mod_seq1)
  
  modification_index_list <- list()
  modification_index_list$pep_seq <- naked_peptide
  modification_index_list$naked_peptide_length <- nchar(naked_peptide)
  
  for ( mod in unique(modification_labels) ) {
    
    ## ## Remove other modifications from the sequence 
    ## that are not being accessed to get a more accurate position index
    
    mods_to_remove_from_sequence <- modification_labels[ !(modification_labels %in% mod) ]
    
    ## In case of 2 PTMs allowed
    if( length(mods_to_remove_from_sequence)==1 ){
      
      current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence )), "", mod_seq1 )
    }else if( length(mods_to_remove_from_sequence)==2 ){
      mid_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence[1] )), "", mod_seq1 )
      current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence[2] )), "", mid_mod_sequence )
    }else{current_mod_sequence <-mod_seq1}
    
    ## Replace current Modification string with an identifier to get amino acid index
    current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mod )), "@", current_mod_sequence )
    
    ## Get positions of current modification based on identifier
    pos_mod <- gregexpr( '@', current_mod_sequence )
    
    
    if (mod == "(p)"){
      if ( as.numeric(unlist((pos_mod)))[1] == 1 ){
        pos_mod <- pos_mod + 1
        
      }else if (length(as.numeric(unlist((pos_mod))))==1){
        pos_mod <-as.numeric(unlist((pos_mod))) 
        
      }else if (length(as.numeric(unlist((pos_mod))))==2){
        
        pos_mod <-as.numeric(unlist((pos_mod))) +1
        pos_mod[2] <- pos_mod[2] +2
        
        
      }else if(length(as.numeric(unlist((pos_mod))))==3){
        
        pos_mod <- as.numeric(unlist((pos_mod))) +1
        
        pos_mod[2] <- pos_mod[2] +2
        pos_mod[3] <- pos_mod[3] +3
      }
    }else if (mod == "(ox)"){
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
    }
    ## Get actual amino acid index
    
    
     
    
    ## Store position in list
    modification_index_list[[paste0('modification_',mod)]] <- as.numeric(unlist(pos_mod))
    
  }
  return( modification_index_list )    
  
}
