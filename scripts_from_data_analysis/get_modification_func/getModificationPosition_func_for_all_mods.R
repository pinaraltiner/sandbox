#
mod_seq <- "AVGM(UniMod:35)PSPVS(UniMod:21)PK"
getModificationPosition_ <- function( mod_seq, character_index=F ){ 
  DEBUG = F
  if ( DEBUG ){
    
    mod_seq <- "EGHAQNPMEPSVPQLS(UniMod:21)LMDVK"
    mod_seq <- "EGHAQNPMEPS(UniMod:21)VPQLS(UniMod:21)LM(UniMod:35)DVK"
    mod_seq <- "ANS(UniMod:21)SPTTNIDHLK(Label:13C(6)15N(2))"
  }
  
  ## Get rid of heavy labels. @TODO fix this later, make more flexible
  mod_seq <- gsub('\\(Label:13C\\(\\d+\\)15N\\(\\d+\\)\\)', '', mod_seq)
  
  if (character_index==F){
    modification_labels <- base::regmatches(mod_seq, gregexpr("\\(.*?\\)", mod_seq))[[1]]
    naked_peptide <- gsub( paste(gsub('\\)','\\\\)',gsub('\\(','\\\\(',modification_labels)), collapse = '|'), '', mod_seq )
    modification_index_list <- list()
    modification_index_list$naked_peptide_length <- nchar(naked_peptide)
    modification_index_list$pep_seq <- naked_peptide
    for ( mod in unique(modification_labels) ) {
      
      ## Remove other modifications from the sequence that are not being accessed to get a more accurate position index
      mods_to_remove_from_sequence <- modification_labels[ !(modification_labels %in%  mod) ]
      if ( length(mods_to_remove_from_sequence)>0 ){
        ## Remove the other modifications not being assessed.
        current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence )), "", mod_seq )
      } else {
        current_mod_sequence <- mod_seq
      }
      ## Replace current Modification string with an identifier to get amino acid index
      current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(",  mod )), "@", current_mod_sequence )
      ## Get positions of current modification based on identifier
      
      
      if(mod == "(Unimod:1)"){
        current_mod_sequence <- gsub("[\t\n]", "",current_mod_sequence)
        pos_mod <- gregexpr( '@', current_mod_sequence )
        
      }else{
        pos_mod <- gregexpr( '@', current_mod_sequence )
      }
      ## Get actual amino acid index
      
      #pos_mod <- pos_mod[[1]] -1
      if(pos_mod[[1]][1] ==1){
        pos_mod[[1]][1] <- pos_mod[[1]][1]
        
      }else{
        
        pos_mod <- pos_mod[[1]] -1
      }
      
      if (length(pos_mod)==2){
        
        pos_mod[2] <- pos_mod[2] -1
        
      }else if(length(pos_mod)==3){
        pos_mod[2] <- pos_mod[2] -1
        pos_mod[3] <- pos_mod[3] -2
      }
      
      ## Store position in list
      modification_index_list[[paste0('modification_',mod)]] <- pos_mod
    }
    return( modification_index_list )    
    
  } else {
    pos_mod <- gregexpr('\\(UniMod:21\\)|\\(Phospho\\)', mod_seq)
    pos_mod <- as.numeric(pos_mod[[1]])-1
  }
  
  ## If there is no modification, pos_mod is negative, then return 0 for pos_mod
  if ( pos_mod < 0 ){
    pos_mod <- 0
  }
  return(pos_mod) 
  
}
