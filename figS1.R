#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure S1

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("gridExtra")

source("setup_env.R")
set.seed(778)


##### LOAD DATA #####
source("load_data.R")


##### SETTINGS #####
save.plots <- TRUE


##### FIG S1 FUNCTIONS #####

##### PROCESS DATASETS #####
clust.df <- list()
for(s in species){
  # Remove non-orthologous DNMs
  clust.df[[s]] <- dnm.df[[s]][, .(CLUST_SIZE = CLUST_SIZE[1], 
                                   TRANSMIT   = all((TRANSMIT_F2A %in% 1) | (TRANSMIT_F2B %in% 1)),
                                   NO_F2      = is.na(F2A[1]) & is.na(F2B[1]),
                                   BLOCK      = BLOCK[1]), by=CLUST_ID]
}

##### PLOT FREQUENCY OF CLUSTERS #####
gg.obj <- list()
for(s in species){
  plt.df <- clust.df[[s]]
  brks <- 1:max(plt.df$CLUST_SIZE)
  gg.obj[[s]] <- ggplot(data=plt.df) + geom_histogram(aes(x=CLUST_SIZE), binwidth=1, fill=sp.color.vec[s], size=1, color="#333333") 
  gg.obj[[s]] <- gg.obj[[s]] + scale_y_log10(breaks=c(1,10,100,1000), limits=c(1,1000)) + annotation_logticks(sides="l") + ylab("Count")
  gg.obj[[s]] <- gg.obj[[s]] + scale_x_continuous(breaks=brks, labels=c("singlet", as.character(brks)[-1]), limits=range(brks) + c(-.5,.5), name="Cluster Size")
  gg.obj[[s]] <- gg.obj[[s]] + ggtitle(sprintf("%s\n", species.nms[s]))
  gg.obj[[s]] <- gg.obj[[s]] + bw.theme + axis.theme + theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank())
}
gg.multi <- grid.arrange(gg.obj[["hum"]], gg.obj[["bab"]], ncol=2, widths=c(6, 8/3))
print(gg.multi)

out.fn <- JoinPath(res.dir, "figS1A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.multi, filename=out.fn, width=6+8/3, height=4, units="in", useDingbats=FALSE)
}

##### PLOT TRANSMISSION OF CLUSTERS #####
gg.obj <- list()
for(s in species){
  cur.df <- clust.df[[s]][!(NO_F2),]
  
  # Expected probability of transmitting a DNM (mean across trios)
  exp.xmit <- mean(tri.df[s,][!duplicated(F1) & !is.na(F2), ComputeQ(Pat, Xmit_Prob), by=F1]$V1)
  
  brks <- 1:max(cur.df$CLUST_SIZE)
  x <- sort(unique(cur.df$CLUST_SIZE))
  y <- c(); upr <- c(); lwr <- c()
  
  for(i in x){
    num <- cur.df[(CLUST_SIZE == i), sum(TRANSMIT)]
    den <- cur.df[(CLUST_SIZE == i), .N]
    y <- c(y, num / den)
    if(num > 0){
      ci <- cur.df[(CLUST_SIZE == i), BootCIByGroup(BLOCK, TRANSMIT)][,"TRUE"]
    } else {
      ci <- c(0,0)
    }
    lwr <- c(lwr, ci[1])
    upr <- c(upr, ci[2])
  }
  plt.df <- data.frame("clust.size"=x, "prop.transmitted"=y, "lwr"=lwr, "upr"=upr)
  gg.obj[[s]] <- ggplot() + geom_pointrange(data=plt.df, mapping=aes(x=clust.size, y=prop.transmitted, ymin=lwr, ymax=upr), color=sp.color.vec[s], shape=21, size=1, fatten=2, fill="white")
  gg.obj[[s]] <- gg.obj[[s]] + geom_hline(yintercept=exp.xmit, color="red", alpha=0.5, size=1)
  gg.obj[[s]] <- gg.obj[[s]] + scale_x_continuous(breaks=brks, labels=c("singlet", as.character(brks)[-1]), name="Cluster Size", limits=range(brks)+c(-.5,.5)) + ylab("Proportion Transmitted")
  gg.obj[[s]] <- gg.obj[[s]] + ylim(0, 0.6)
  gg.obj[[s]] <- gg.obj[[s]] + ggtitle(species.nms[s])
  gg.obj[[s]] <- gg.obj[[s]] + bw.theme + axis.theme + theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
}
gg.multi <- grid.arrange(gg.obj[["hum"]], gg.obj[["bab"]], ncol=2, widths=c(6, 8/3))

print(gg.multi)
out.fn <- JoinPath(res.dir, "figS1B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.multi, filename=out.fn, width=6+8/3, height=3.8, units="in", useDingbats=FALSE)
}
