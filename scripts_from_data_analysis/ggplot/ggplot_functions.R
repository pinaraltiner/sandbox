#
library(ggplot2)
#library(dplyr)
library(gghalves)

##Distribution of both mean abundance 
#                                    and Exp. Ratios 
#                                                   of every sample with two pools
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
    scale_fill_brewer(palette = 2,direction=-1) +
    scale_linetype_manual(values=c("dashed", "dotted")) +
    theme_minimal() +
    theme(legend.text = element_text(size=25), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25),
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
#         axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),
#         plot.title = element_text(size=30),
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
    theme(legend.text = element_text(size=25), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)
    ) +   stat_boxplot(geom = "errorbar") + 
    ggtitle(header) +
    labs(x=x_lab,y=y_lab,fill = fill_lab,subtitle = subtitle_txt) +
    scale_fill_brewer(palette="Set1")
  
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
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)
    ) + scale_fill_brewer(palette="Set1") +
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
    theme(legend.text = element_text(size=25),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)
    ) + scale_fill_brewer(palette="Set1") +
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
                                subtitle_txt){


ggplot(data_set, aes(x=x_df, fill=fill_df)) + geom_bar(position = "dodge") +
  scale_fill_brewer(palette = 2,direction=-1) +
  theme_minimal() +
  theme(legend.text = element_text(size=25), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=30),
        plot.subtitle = element_text(size = 20),
        plot.caption = element_text(size = 20),
        legend.title=element_text(size=25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=25)) +
  ggtitle(header) + 
  scale_y_continuous(breaks = seq(from=0, to=ymax,by=500)) +
  geom_text(aes(label=after_stat(count)),stat = "count", position=position_dodge(width=0.9), vjust=-0.25) +
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
                                         y = log10(y_df),
                                         fill=factor(fill_df))) +
    ggdist::stat_halfeye( adjust = 0.5,
                          #justification = -2,
                          .width = 0,
                          point_colour=NA,
                          alpha=0.8) +
    geom_boxplot(width = .12, 
                 outlier.colour = NA,
                 alpha=0.5) +
    ggdist::stat_dots(side="left",
                      justification= 1.1,
                      binwidth =.025) + 
    theme_light() + 
    theme(legend.text = element_text(size=25), #plot.margin=unit(c(-0.5,1,1,1), "cm"),
          axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),
          plot.title = element_text(size=30),
          plot.subtitle = element_text(size = 20),
          plot.caption = element_text(size = 20),
          legend.title=element_text(size=25),
          axis.text=element_text(size=25),
          axis.title=element_text(size=25)) + 
    ggtitle(header) + 
    labs(x=x_lab,y=y_lab, fill= fill_lab, caption = caption_lab,subtitle = subtitle_txt) +
    scale_fill_brewer(palette="Set1")
  
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
      scale_y_continuous(limits = c(round(min(y_df),2), round(max(y_df),2)), breaks = seq(round(min(y_df),2), round(max(y_df),2), by = 1)) +
      scale_x_continuous(limits = c(min(x_df), max(x_df)),breaks = seq(min(x_df), max(x_df), by = 1)) +
      scale_color_brewer(palette = "Set1") +
      #scale_y_continuous(breaks = seq(0, max(-log10(volcano_final1$pvalues_value)), length.out = 21)) +
      theme_bw() +
      theme(legend.text = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 30),
          legend.title = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    #expand_limits(x = 0, y = 0) +
    #geom_vline(data = actual_ratio, aes(xintercept = actual_ratio$X.1....log2.c.2..10..20..100..., size = 1, show.legend = FALSE)) + #color=c("#CC79A7","#E69F00","#56B4E9","#009E73")
    #geom_hline(data = log10_p_thresholds, aes(yintercept = log10_p_thresholds$X.log10.p_thresholds.),color=c("#CC79A7","#E69F00","#56B4E9","#009E73"), size = 1, linetype = 2, show.legend = FALSE)+ 
    labs(title = header, color = color_lab, subtitle = subtitle_txt)}
  
}
  
  
  
  
  
  