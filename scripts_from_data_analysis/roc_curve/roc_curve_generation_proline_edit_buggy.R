#' Compute ROC curve based on the pvalue from the supplied dataset. pvalue and flag columns must have been filled before
#' calling this function.  
#'
#' @param df The dataset
#' @param flag The flag value identifying expected differentially expressed proteins  
#' @param expected The theoritically expected number of differential proteins
#'
#' @return
#' @export
#'
#'

compute_roc_curve = function(df, flag, expected) {
  data <- df %>% 
    rename(pvalue=contains("Value")) %>%
    #arrange_at(vars(starts_with("p_val"))) %>%
    #rename(pvalue = starts_with("p_val")) %>%
    rename(flag=contains("Pool")) %>% drop_na(flag)
  result = data.frame(fdp=double(nrow(data)-1),tpr=double(nrow(data)-1),pvalue = double(nrow(data)-1), stringsAsFactors = F)
  
  
  tp = nrow(data[data$flag == flag,])
  fp = nrow(data[data$flag != flag,])
  for (k in nrow(data):2) {
    if (data$flag[k] == flag) {
      tp = tp - 1
      #print(tp)
    } else if (data$flag[k] != flag) {
      fp = fp - 1
      #print(fp)
    }
    
    fdp = fp#(fp/(fp+tp))*100
    tpr = tp#(tp/expected)*100 
    
    result$fdp[nrow(data) - (k-1)] = fdp
    result$tpr[nrow(data) - (k-1)] = tpr
    
    result$pvalue[nrow(data) - (k-1)] = data$pvalue[k]
  }
  
  lastValue = result$fdp[1]
  for (i in 2:(nrow(result))) {
    if (!is.na(result$fdp[i]) & (result$fdp[i] > lastValue)) { 
      result$fdp[i] = NA
    } else {
      lastValue = result$fdp[i]
    }  
  }
  result = result[!is.na(result$fdp),]
  result = result[order(result$tpr, result$fdp, decreasing = TRUE),]
  return(result)
}



# 
# compute_roc_curve = function(df, flag, expected) {
#   data <- df %>% 
#     rename(pvalue=contains("Value")) %>%
#     #arrange_at(vars(starts_with("p_val"))) %>%
#     #rename(pvalue = starts_with("p_val")) %>%
#     rename(flag=contains("Pool")) %>% drop_na(flag)
#   
#   
#   positives <- data %>% filter(pvalue < 0.05)
#   negatives <- data %>% filter(pvalue > 0.05)
#   
#   result = data.frame(fdp=double(nrow(positives)-1),tpr=double(nrow(positives)-1),pvalue = double(nrow(positives)-1), stringsAsFactors = F)
#   tp = nrow(positives[positives$flag == flag,])
#   fp = nrow(positives[positives$flag != flag,])
#   
#   for (k in nrow(positives):2) {
#     if (positives$flag[k] == flag) {
#       tp = tp - 1
#       #print(tp)
#     } else if (positives$flag[k] != flag) {
#       fp = fp - 1
#       #print(fp)
#     }
#     
#     fdp = (fp/(fp+tp))*100
#     tpr = (tp/expected)*100 
#     
#     result$fdp[nrow(positives) - (k-1)] = fdp
#     result$tpr[nrow(positives) - (k-1)] = tpr
#     
#     result$pvalue[nrow(positives) - (k-1)] = positives$pvalue[k]
#   }
#   
#   lastValue = result$fdp[1]
#   for (i in 2:(nrow(result))) {
#     if (!is.na(result$fdp[i]) & (result$fdp[i] > lastValue)) { 
#       result$fdp[i] = NA
#     } else {
#       lastValue = result$fdp[i]
#     }  
#   }
#   result = result[!is.na(result$fdp),]
#   result = result[order(result$tpr, result$fdp, decreasing = TRUE),]
#   return(result)
# }
# 
