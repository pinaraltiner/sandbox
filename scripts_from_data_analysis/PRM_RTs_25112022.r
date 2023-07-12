library(ggplot2)
library(hrbrthemes)
prm_rts <- read.table("RT Errors.txt", sep = "\t", header = T)
colnames(prm_rts) <- c("before RT calibration","after RT calibration")
d <- density(prm_rts$DELTA_PRM_DDA)
d1 <- density(prm_rts$RT_ERROR)
plot(d, main="",col="red", ylim=c(0,0.06),xlim=c(-100,100))
par(new=TRUE)
plot(d1, main="",col="blue",ylim=c(0,0.06),xlim=c(-100,100))
#polygon(d, col="red", border="blue")
"before RT calibration"
"after RT calibration"

legend(x = "topleft", legend=c("before RT calibration",
                               "after RT calibration"), 
       fill = c("red","blue")
)

library(reshape)
melt_df <- melt(prm_rts)

p2 <- ggplot(data=melt_df, aes(x=value,group=variable, fill=variable)) +
    geom_density(adjust=1.5, alpha=.4, size=1) +  scale_x_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
    theme_ipsum() +
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
          legend.text = element_text(size=15),
          plot.margin=unit(c(1,1,-0.5,1), "cm"),
          plot.title = element_text(size=20), #legend.position = "none",
    ) + labs(x="Retention Time error",y="Density") +
    guides(fill=guide_legend(title="")) +
    ggtitle("Distribution of Retention Time (RT) errors")
p2
ggsave("dist_of_rt_before_after_calib.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

library(viridis)
library(patchwork)

p <- ggplot(melt_df, aes(x="", y=value, fill=variable)) +
    geom_boxplot(alpha=0.7) + coord_flip() +
    #stat_summary(funy =mean, geom="point", shape=20, size=14, color="red", fill="red") +
    theme_ipsum() + scale_y_continuous(limits = c(-60, 60),breaks = seq(-60, 60, by = 20)) +
    theme(plot.margin=unit(c(-0.5,1,1,1), "cm"),legend.text = element_text(size=15),
      axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
      plot.title = element_text(size=20),
    ) +  guides(fill=guide_legend(title="")) + stat_boxplot(geom = "errorbar") +
    ggtitle("") +
    labs(x="",y="Retention Time error (s)")
    #scale_fill_brewer(palette="Set1")
p2/p
ggsave("new_closer_denisty_boxplot_horizontal_two_legend_dist_of_rt_before_after_calib.tiff", units="in", width=9, height=6, dpi=300, compression = 'lzw') #w=9.5, h=5



