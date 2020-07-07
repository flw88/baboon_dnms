#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure 3

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("digest")
library("ggrepel")

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


##### CALCULATE ALPHA FROM GLM #####
# Load regression results
glm.fn <- JoinPath(res.dir, "glm.ortho.RData")
StopNoGLMData(glm.fn)

load(glm.fn) # Loads save.list
glm.obj <- save.list$glm.obj
glm.sim <- save.list$glm.sim

# Calculate alpha
alph.glm <- c()
alph.glm.ci <- list()
for(s in species){
  
  mut.counts <- c()
  sim.counts <- c()
  for(p in parents){
    mut.counts[p] <- predict(glm.obj[[s]][[p]], newdata=data.frame("age"=avg.gen[[s]][p]))
    sim.counts <- cbind(sim.counts, colSums(t(coef(glm.sim[[s]][[p]])) * c(1, avg.gen[[s]][p])))
  }
  colnames(sim.counts) <- parents
  
  alph.glm[s] <- mut.counts["Pat"] / sum(mut.counts)
  alph.glm.ci[[s]] <- quantile(sim.counts[,"Pat"] / rowSums(sim.counts), probs=c(.025,.975))
  names(alph.glm.ci[[s]]) <- c("lwr","upr")
}
  


##### CALCULATE ALPHA FROM TRANSMITTED DNMS #####
CalcAlpha <- function(df.list, mask.list=NULL){
  phase.xmit <- matrix(0, nrow=2, ncol=2, dimnames= list(species, parents))
  alph  <- data.table("sp"=species)
  setkey(alph, sp)
  
  for(s in names(df.list)){
    cur.df <- df.list[[s]]
    if(!is.null(mask.list)){
      cur.df <- cur.df[mask.list[[s]],]
    }
    
    xmit.mask <- (cur.df$TRANSMIT_F2A %in% 1) | (cur.df$TRANSMIT_F2B %in% 1)
    phase.vec <- cur.df[xmit.mask, PHASE_CONS]
    grp.vec   <- cur.df[xmit.mask, BLOCK]
    
    m <- phase.vec %in% p.chars
    phase.vec <- phase.vec[m]
    grp.vec   <- grp.vec[m]
    
    # Calculate confidence intervals
    for(p in parents){
      phase.xmit[s,p] <- sum(phase.vec == p.chars[p], na.rm=T)
    }
    
    
    alph[s, c("lwr", "upr") := as.list(BootCIByGroup(grp.vec, phase.vec)[,"P"]) ]
  }
  
  # Calculate alpha
  alph[, "a" := phase.xmit[sp,1] / phase.xmit[sp,2]]
  alph[, "prop" := a / (1+a) ]
  
  return(list("alph"=alph, "phase.xmit"=phase.xmit))
}

mask.list <- list()
for(s in species){
  mask.list[[s]] <- dnm.df[[s]][, TYPE != "CpG>TpG"]
}

tmp.all    <- CalcAlpha(dnm.df)
tmp.noncpg <- CalcAlpha(dnm.df, mask.list) # Without CpG Ti's

phase.xmit <- list("all"    = tmp.all$phase.xmit,
                   "noncpg" = tmp.noncpg$phase.xmit)
alph <- list("all"    = tmp.all$alph,
             "noncpg" = tmp.noncpg$alph)



# Calculate p-value
test.res <- list()
for(nm in names(phase.xmit)){
  test.res[[nm]] <- chisq.test(phase.xmit[[nm]])
}

if(save.plots){
  alph.stats.fn <- JoinPath(res.dir, "report.alpha.txt")
  cat("Saving", alph.stats.fn, "\n")
} else {
  alph.stats.fn <- "/dev/stdout"
}

cat("========== ALPHA STATISTICS ==========\n", file=alph.stats.fn, append=FALSE)
# Print phase.xmit to stdout
for(nm in names(alph)){
  cat("### Using", nm, "dnms ###\n", file=alph.stats.fn, append=TRUE)
  cat("\nPhase Table:\n", file=alph.stats.fn, append=TRUE)
  sink(alph.stats.fn, append=TRUE)
  print(phase.xmit[[nm]])
  cat("\nAlpha:\n", file=alph.stats.fn, append=TRUE)
  print(alph[[nm]])
  cat(sprintf("P-value for sp. diff = %0.3f \n", test.res[[nm]]$p.value), file=alph.stats.fn, append=TRUE)
  sink()
  cat("\n", file=alph.stats.fn, append=TRUE)
}
  

##### PLOT ALPHA #####

# Plot proportion paternal
plt.df <- data.frame("x"=seq_along(species)-1/8,
                     "y"=phase.xmit$all[,"Pat"]/rowSums(phase.xmit$all),
                     "species"=species, "method"="transmit")
plt.df <- cbind(plt.df, alph$all[species, .(lwr, upr)])
for(i in seq_along(species)){
  s <- species[i]
  tmp <- data.frame("x"=i+1/8, "y"=alph.glm[s], 
                    "species"=s, "method"="regression",
                    "lwr"=alph.glm.ci[[s]]["lwr"],
                    "upr"=alph.glm.ci[[s]]["upr"])
  plt.df <- rbind(plt.df, tmp)
}

gg.obj <- ggplot(data=plt.df) + geom_pointrange(aes(x=x, y=y, ymin=lwr, ymax=upr, color=species, shape=method), size=1, fatten=2, fill="white")
gg.obj <- gg.obj + scale_shape_manual(values=c(21,22))
gg.obj <- gg.obj + scale_x_continuous(name="", breaks=seq_along(species), labels=species.nms[species], limits=c(0.5, length(species)+0.5))
gg.obj <- gg.obj + scale_color_manual(values=sp.color.vec[levels(plt.df$species)])
gg.obj <- gg.obj + ylim(0,1) + ylab(bquote(Proportion ~ Paternal ~ (alpha))) + guides(fill=F)
gg.obj <- gg.obj + ggtitle("Male mutation bias estimated from\ntransmitted, phased DNMs")
gg.obj <- gg.obj + bw.theme + axis.theme + xgrid.off.theme
print(gg.obj)

out.fn <- JoinPath(res.dir, "fig3A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=4, height=4, units="in", useDingbats=FALSE)
}


##### PLOT ALPHA VS. GEN TIME FOR VARIETY OF SPECIES #####

alpha.tab["Homo sapiens","Ratio"] <- alph$all["hum",a]
alpha.tab["Papio anubis","Ratio"] <- alph$all["bab",a]

alpha.tab[, y := Ratio / (1+Ratio)]
alpha.tab[, Color := "grey15"]
alpha.tab["Homo sapiens", Color := sp.color.vec["hum"]]
alpha.tab["Papio anubis", Color := sp.color.vec["bab"]]
alpha.tab[, Color := as.factor(Color)]

### Perform a linear regression
RegrAlpha <- function(alpha.tab, rm.sp = c(), out.fn){
  lm.tab <- alpha.tab[!(Species %in% rm.sp),]
  alpha.lm <- lm(y ~ Mean_Pat_Age, data=lm.tab)
  
  # Print 
  summ <- coef(summary(alpha.lm))
  coef.pt <- coef(alpha.lm)
  coef.ci <- confint(alpha.lm)
  
  sink(out.fn, append=TRUE)
  cat("\n====================================================\n", file=out.fn, append=TRUE)
  cat("Linear regr. of alpha vs. pat. age after removing\n", paste(rm.sp, collapse=",\n"), ":\n", sep="", file=out.fn, append=TRUE)
  print(summ)
  sink()
  cat("\n", file=out.fn, append=TRUE)
  # Additionally print intercept and slope with CI
  
  fmt.str <- "%15s   %0.1e (%0.1e -- %0.1e)   %0.1e\n"
  for(nm in names(coef.pt)){
    cat(sprintf(fmt.str, nm, coef.pt[nm], coef.ci[nm,1], coef.ci[nm,2], summ[nm,4]), file=out.fn, append=TRUE)
  }
  
  return(list("lm.tab"=lm.tab, "alpha.lm"=alpha.lm))
}


rm.sp <- c("Pan troglodytes (Tatsumoto et al.)", "Homo sapiens (Jónsson et al.)")
tmp <- RegrAlpha(alpha.tab, rm.sp, alph.stats.fn)
lm.tab <- tmp$lm.tab
alpha.lm <- tmp$alpha.lm

# Also print what would happen if we removed the other chimpanzee estimate.
tmp <- RegrAlpha(alpha.tab, c("Pan troglodytes (Besenbacher et al.)", rm.sp[2]), alph.stats.fn)


# Build data tables for plot 
tmp <- data.frame("Mean_Pat_Age"=SeqRange(lm.tab[["Mean_Pat_Age"]]))
lm.pred <- predict.lm(alpha.lm, tmp, se.fit=TRUE)
line.df <- data.table("x"=tmp$Mean_Pat_Age, "y"=lm.pred$fit, 
                      "lwr"=lm.pred$fit + lm.pred$se.fit * qnorm(0.025),
                      "upr"=lm.pred$fit + lm.pred$se.fit * qnorm(0.975))


tab.rn <- alpha.tab[, Species]
use.shapes <- rep(21, nrow(alpha.tab)) # Circle
use.shapes[tab.rn == "Pan troglodytes (Tatsumoto et al.)"] <- 23 # Diamond
use.shapes[tab.rn == "Pan troglodytes (Besenbacher et al.)"] <- 22 # Square
use.shapes[tab.rn == "Homo sapiens (Jónsson et al.)"] <- 25 # Triangle
use.shapes <- as.factor(use.shapes)

# Nudge the owl monkey point
nudge.x <- rep(0, nrow(alpha.tab))
nudge.x[tab.rn == "Aotus nancymaae"] <- -3

nudge.y <- rep(0, nrow(alpha.tab))
nudge.y[tab.rn == "Aotus nancymaae"] <- 0


gg.obj <- ggplot(data=alpha.tab)
gg.obj <- gg.obj + geom_line(aes(x=x, y=y), data=line.df, color="grey35", linetype="dashed")
gg.obj <- gg.obj + geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), data=line.df, fill="black", alpha=0.2)

gg.obj <- gg.obj + geom_point(aes(x=Mean_Pat_Age, y=y, color=Color, shape=use.shapes, fill=Color), size=2.1, show.legend=FALSE)
gg.obj <- gg.obj + scale_shape_manual(values=as.integer(levels(use.shapes)))
gg.obj <- gg.obj + geom_text_repel(aes(x=Mean_Pat_Age, y=y, color=Color, label=Common), fontface="plain", 
                                   size=3.5, show.legend=FALSE, force=7, seed=3, nudge_x=nudge.x, nudge_y=nudge.y)

gg.obj <- gg.obj + scale_color_manual(values=levels(alpha.tab$Color))
gg.obj <- gg.obj + scale_fill_manual( values=levels(alpha.tab$Color))
gg.obj <- gg.obj + xlab("Mean Paternal Age at Conception (y)") + ylab("Proportion Paternal DNMs") + xlim(-2,37) + ylim(-0.05,1)
gg.obj <- gg.obj + ggtitle("Alpha and Generation Time Estimates\nfrom Pedigree Studies") + bw.theme + axis.theme
print(gg.obj)

out.fn <- JoinPath(res.dir, "fig3B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=5, height=4, units="in", useDingbats=FALSE)
}
