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

# compute_roc_curve = function(df, flag, expected) {
#   data <- df %>% 
#     rename(pvalue=contains("Value")) %>%
#     #arrange_at(vars(starts_with("p_val"))) %>%
#     #rename(pvalue = starts_with("p_val")) %>%
#     rename(flag=contains("Pool")) %>% drop_na(flag)
#   result = data.frame(fdp=double(nrow(data)-1),tpr=double(nrow(data)-1),pvalue = double(nrow(data)-1), stringsAsFactors = F)
#   
#   
#   tp = nrow(data[data$flag == flag,])
#   fp = nrow(data[data$flag != flag,])
#   for (k in nrow(data):2) {
#     if (data$flag[k] == flag) {
#       tp = tp - 1
#       #print(tp)
#     } else if (data$flag[k] != flag) {
#       fp = fp - 1
#       #print(fp)
#     }
#     
#     fdp = (fp/(fp+tp))*100
#     tpr = (tp/expected)*100 
#     
#     result$fdp[nrow(data) - (k-1)] = fdp
#     result$tpr[nrow(data) - (k-1)] = tpr
#     
#     result$pvalue[nrow(data) - (k-1)] = data$pvalue[k]
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



compute_roc_curve = function(df, flag, expected) {
  data <- df %>%
    rename(pvalue=contains("Value")) %>%
    #arrange_at(vars(starts_with("p_val"))) %>%
    #rename(pvalue = starts_with("p_val")) %>%
    rename(flag_name=contains("pool")) %>% drop_na(flag_name)
  
  thresholds <- thresholds <- c(1E-10,2E-10,4E-10,8E-10,1.6E-09,
                                3.2E-09,6.4E-09,1.28E-08,2.56E-08,
                                5.12E-08,1.024E-07,2.048E-07,
                                4.096E-07,8.192E-07,1.6384E-06,
                                3.2768E-06,6.5536E-06,1.31072E-05,
                                2.62144E-05,5.24288E-05,0.000104858,
                                0.000209715,0.00041943,0.000838861,
                                0.001677722,0.003355443,0.006710886,
                                0.013421773,0.026843546,0.053687091,
                                0.107374182,0.214748365,0.42949673,
                                0.858993459)
    #c(0.060,0.055,0.050,0.045, 0.040, 0.035,0.030,0.025, 0.020,0.015, 0.010)
  result <- NULL
  for (i in 1:length(thresholds)){
    
    tp <- dim(as.data.frame(data %>% filter(pvalue < thresholds[i] & grepl(flag,flag_name))))[1]
    fp <- dim(as.data.frame(data %>% filter(pvalue < thresholds[i] & !grepl(flag,flag_name))))[1]
    
    #tp = positives[positives$flag == flag,])
    #fp = positives[positives$flag != flag,])
    
    fdr = (fp/(fp+tp))*100
    tpr = (tp/expected)*100
    
    result$tpr[i] <- tpr
    result$fdr[i] <- fdr
    
  }
  return(as.data.frame(result))
  
  
# 
#   for (k in nrow(positives):2) {
#     if (positives$flag[k] == flag) {
#       tp = tp - 1
#       #print(tp)
#     } else if (positives$flag[k] != flag) {
#       fp = fp - 1
#       #print(fp)
#     }

  #   fdp = (fp/(fp+tp))*100
  #   tpr = (tp/expected)*100
  # 
  #   result$fdp[nrow(positives) - (k-1)] = fdp
  #   result$tpr[nrow(positives) - (k-1)] = tpr
  # 
  #   result$pvalue[nrow(positives) - (k-1)] = positives$pvalue[k]
  # }

  # lastValue = result$fdp[1]
  # for (i in 2:(nrow(result))) {
  #   if (!is.na(result$fdp[i]) & (result$fdp[i] > lastValue)) {
  #     result$fdp[i] = NA
  #   } else {
  #     lastValue = result$fdp[i]
  #   }
  # }
  # result = result[!is.na(result$fdp),]
  # result = result[order(result$tpr, result$fdp, decreasing = TRUE),]
}

