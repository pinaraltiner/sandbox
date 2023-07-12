#
library('ggplot2')
library('opentimsr')
library("reshape2")
path ="D:/dev/Desktop_copy/PHD/wet_lab_experiments/secondment_sample_prep/data_from_sdu/all_transfer_files/"

folder_id <- c("20230302-1517_TTP_002917_ONJ_TR_Toulouse_E1_M1_M8coli_pool_CE_5420.d",
  "20230302-1428_TTP_002916_ONJ_TR_Toulouse_E1_M1_M8_pool_CE_5419.d",
  "20230301-1913_TTP_002912_ONJ_TR_Toulouse_E1_M1_M8coli_pool_5415.d",
  "20230301-1825_TTP_002911_ONJ_TR_Toulouse_E1_M1_M8_pool_5414.d"
  )

accept_Bruker_EULA_and_on_Windows_or_Linux = TRUE
  
  if(accept_Bruker_EULA_and_on_Windows_or_Linux){
    
    folder_to_stode_priopriatary_code = "D:/dev/Desktop_copy/PHD/wet_lab_experiments/secondment_sample_prep/data_from_sdu/"
    path_to_bruker_dll = download_bruker_proprietary_code(folder_to_stode_priopriatary_code)
    setup_bruker_so(path_to_bruker_dll)
    all_columns = c('frame','scan','tof','intensity','mz','inv_ion_mobility','retention_time')
    
  } else {
    
    all_columns = c('frame','scan','tof','intensity','retention_time')
  }

for (i in 1:length(folder_id)){
  
  path_tmp <- paste0(path,folder_id[i])
  
  
  D = OpenTIMS(path_tmp)
  assign(paste0("result_mz",i), NULL)
  assign(paste0("result_int",i), NULL)
  all_final_comb <- NULL
  assign(paste0("all_prec_count",i), NULL)
  # All MS1 frames, but one at a time:
  
  for(fr in MS1(D)){
    # number obtained precursor for each frame
    #tmp <- dim(query(D, fr, columns=all_columns[c(2,5)]))
    
    # Median of number of identified precursor for each frame
    tmp2 <- median(as.numeric(unlist(query(D, fr, columns=all_columns[c(4,5)])[2])))
    tmp3 <- median(as.numeric(unlist(query(D, fr, columns=all_columns[c(4,5)])[1])))
    assign(paste0("result_mz",i),rbind(get(paste0("result_mz",i)),tmp2))
    assign(paste0("result_int",i),rbind(get(paste0("result_int",i)),tmp3))
    #assign(paste0("all_prec_count",i),rbind(get(paste0("all_prec_count",i)),dim(tmp)))
  

  }
  rm(D)
}

#  SAME AS all_prec_count_1,2,3,4 
num_ms1 <- t(cbind(dim(result_mz1),dim(result_mz2),dim(result_mz3),dim(result_mz4)))
rownames(num_ms1) <- folder_id  


ggplot(data = as.data.frame(num_ms1), aes(y=V1, x=1:4, fill=rownames(num_ms1))) +
  geom_bar(stat = "identity") + theme_bw() +
  theme_bw() + coord_flip() + 
  theme(legend.text = element_text(size=15), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size=20),
        legend.title=element_text(size=15),
        axis.text=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(size=15)
  ) + scale_y_continuous(breaks = seq(0,max(num_ms1)+100,700)) +
  labs(x="number of precursor ions",y="methods")


for (i in 1:length(folder_id)){
  # Merging mz and int values 
  assign(paste0("final_result",i),cbind(get(paste0("result_mz",i),get(paste0("result_int",i))),folder_id[i]))
  assign(rownames(get(paste0("final_result",i))),NULL)
  assign(get(paste0("all_final_comb")),rbind(all_final_comb,get(paste0("final_result",i)))) 
  
}  



# REMOVE THEM BEFORE CLOSING -> THIS IS FOR ONLY ONCE
final_result1 <- cbind(result_int1,result_mz1,folder_id[1])
final_result2 <- cbind(result_int1,result_mz1,folder_id[2])
final_result3 <- cbind(result_int1,result_mz1,folder_id[3])
final_result4 <- cbind(result_int1,result_mz1,folder_id[4])


# REMOVE THEM BEFORE CLOSING -> THIS IS FOR ONLY ONCE
rownames(final_result1) <- NULL
rownames(final_result2) <- NULL
rownames(final_result3) <- NULL
rownames(final_result4) <- NULL

# REMOVE THEM BEFORE CLOSING -> THIS IS FOR ONLY ONCE
all_final_comb <- rbind(final_result1,final_result2,final_result3,final_result4)


colnames(all_final_comb) <- c("intensity","mz", "file_name")

melt_all_final_comb <- melt(all_final_comb,id.vars=file_name)

hist(final_result1[,2], breaks = 100,xlim = c(200, max(final_result1[,2])))

library(dplyr)

as.data.frame(all_final_comb) %>% group_by(file_name) %>% ggplot(aes(intensity,
                                         color=file_name))+
         geom_histogram(binwidth=10, boundary=-0.5, stat = "count") + 
  scale_x_continuous(breaks=break_leng_mz) + theme_bw()



as.data.frame(all_final_comb) %>% group_by(file_name) %>% ggplot(aes(mz,
                                                                     color=file_name))+
  geom_density(stat="count") + 
  scale_x_continuous(breaks=break_leng_mz) + theme_bw()
  
break_leng_mz <- seq(from=200, to =max(all_final_comb[,2]),by =20)
break_leng_int <- seq(from=0, to=max(as.numeric(all_final_comb[,1])),by=10)

ggplot(as.data.frame(all_final_comb),aes(x =intensity,color=file_name))+geom_density()+ 
  scale_x_discrete(breaks = break_leng_int) + theme_bw()
  




for(fr in MS1(D)){
  tmp <- query(D, fr, columns=all_columns[c(2,5)])
  result <- tmp
}

print(length(D))
pprint = function(x,...){ print(head(x,...)); print(tail(x,...)) }

# Get a data,frame with data from frames 1, 5, and 67.
tmp <- pprint(query(D, frames=c(1,50,100), columns=all_columns))
colnames(tmp)

