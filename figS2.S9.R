#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figures S2 and S9

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("gridExtra")
source("setup_env.R")

##### LOAD DATA #####
source("load_data.R")


##### SETTINGS #####
save.plots <- TRUE


##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]

  # Remove non-orthologous DNMs
  dnm.df[[s]] <- dnm.df[[s]][ORTHO == 1]
}


##### COMPUTE FDR #####
fdr.ones <- list()
fdr.twos <- list()
for(s in species){
  fdr.ones[[s]] <- ComputeSingletFDR(tri.df[s], dnm.df[[s]], ci=0.95)
}

# Estimate FDR for baboon F1s without an F2
s <- "bab"
fdr.lm <- RegressFDR(tri.df[s], fdr.ones[[s]])

##### PLOT FNR RESULTS #####
BootResampleMut <- function(x, n, n.sim=1000, ci=c(0.025, 0.975)){
  prop <- x / n
  sim.vals <- replicate(n.sim, expr=sum(sample(c(0,1), size=n, prob=c(1-prop, prop), replace=TRUE)) / n)
  return(quantile(sim.vals, probs=ci))
}

GGplotErrorRates <- function(plt.df, y.lim=c(0,1), use.color="black"){
  out.obj <- ggplot(plt.df)
  out.obj <- out.obj + geom_pointrange(aes(x=x, y=est, ymin=lwr, ymax=upr), 
                                       shape=21, size=1, fatten=2, fill="white", color=use.color)
  out.obj <- out.obj + scale_x_continuous(name="", breaks=plt.df[,x], labels=plt.df[,F1]) 
  out.obj <- out.obj + ylim(y.lim[1], y.lim[2]) + bw.theme + axis.theme + xgrid.off.theme
  return(out.obj)
}

n.fnr.sim <- 2e4
gg.obj <- list()
for(s in species){
  plt.df <- tri.df[s, .(est=FNR[1]), by=F1]
  cur.ci <- sapply(plt.df$est, function(x) BootResampleMut(x*n.fnr.sim, n=n.fnr.sim))
  plt.df[, lwr := cur.ci[1,]]
  plt.df[, upr := cur.ci[2,]]
  plt.df[, x := 1:.N]
  
  gg.obj[[s]] <- GGplotErrorRates(plt.df, y.lim=c(0, 0.2), use.color=sp.color.vec[s])
  gg.obj[[s]] <- gg.obj[[s]] + ylab("FNR")
  gg.obj[[s]] <- gg.obj[[s]] + ggtitle(paste0(species.nms[s], " False Negative Rates estimated\npedigree simulated DNMs"))
}

gg.multi <- grid.arrange(gg.obj[["hum"]], gg.obj[["bab"]], nrow=2)

out.fn <- JoinPath(res.dir, "figS2A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.multi, filename=out.fn, width=6.7, height=5, units="in", useDingbats=FALSE)
}

##### PLOT FDR RESULTS #####
gg.obj <- list()
for(s in species){
  # Get the species into the right order
  f1.ordered <- unique(tri.df[s,F1])
  plt.df <- data.table::copy(fdr.ones[[s]][f1.ordered][!is.na(lwr)])
  plt.df[, x := 1:.N]
  
  gg.obj[[s]] <- GGplotErrorRates(plt.df, y.lim=c(0,1), use.color=sp.color.vec[s])
  gg.obj[[s]] <- gg.obj[[s]] + ylab("FDR")
  gg.obj[[s]] <- gg.obj[[s]] + ggtitle(paste0(species.nms[s], " False Discovery Rates estimated\nusing transmission data"))
}

gg.multi <- grid.arrange(gg.obj[["hum"]], gg.obj[["bab"]], nrow=2)

out.fn <- JoinPath(res.dir, "figS2B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.multi, filename=out.fn, width=6.5, height=5, units="in", useDingbats=FALSE)
}


##### PLOT FDR REGRESSION #####
cur.df <- fdr.ones[["bab"]]
line.df <- data.table( "depth"= SeqRange(cur.df$depth, length.out=100) )
fit.obj <- predict(fdr.lm, newdata=line.df, se.fit=TRUE)
line.df[, est := fit.obj$fit]
line.df[, lwr := (qnorm(0.025) * fit.obj$se.fit) + est]
line.df[, upr := (qnorm(0.975) * fit.obj$se.fit) + est]

# Regression layers
gg.obj <- ggplot() + geom_line(aes(x=depth, y=est), data=line.df)
gg.obj <- gg.obj + geom_ribbon(aes(x=depth, ymin=lwr, max=upr), data=line.df, 
							   fill="black", alpha=0.2)

# Points layers
gg.obj <- gg.obj + geom_point(aes(x=depth, y=est), data=cur.df[!(pred)],
							  fill=sp.color.vec["bab"], size=3, shape=21, stroke=1)
gg.obj <- gg.obj + geom_point(aes(x=depth, y=est), data=cur.df[(pred)], 
						      fill=NA, size=3, shape=3, stroke=1, color="#3B14AF")
gg.obj <- gg.obj + ylim(-0.3,1.0) + xlab("F1 Mean Depth of Coverage") + 
		  ylab("False Discovery Rate") + bw.theme + axis.theme +
		  ggtitle("Linear Regression of Singlet FDR vs F1 Depth\n")
print(gg.obj)

out.fn <- JoinPath(res.dir, "figS9.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=6, height=4, units="in", useDingbats=FALSE)
}
