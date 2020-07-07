#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plot in Figure S5A

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("arm")

source("setup_env.R")


##### LOAD DATA #####
source("load_data.R")


##### CONSTANTS #####
regr.types <- c("ortho", "all")


##### SETTINGS #####
save.plots <- TRUE


##### LOAD REGRESSION RESULTS #####
glm.obj <- list(); glm.sim <- list()
for(nm in regr.types){
  glm.fn <- JoinPath(res.dir, paste0("glm.", nm, ".RData"))
  StopNoGLMData(glm.fn)
  
  load(glm.fn) # Loads save.list
  glm.obj[[nm]] <- save.list$glm.obj
  glm.sim[[nm]] <- save.list$glm.sim
}


##### PLOT MUT PER GEN WITH CI #####
CalcRateSum <- function(obj.list, age.list, denom){ 
  # Takes in list of "glm" or "sim" objects and ages (names=c("Pat", "Mat"))
  mu.g <- 0
  for(p in names(obj.list)){
    mu.g <- mu.g + CalcRateAtAge(coef(obj.list[[p]]), age.list[[p]], denom)
  }
  
  return(unname(mu.g)) # Return mutation rates for sampled coefficients
}

ci <- 0.95; probs <- c(1-ci, 1+ci) / 2
plt.df <- data.frame()
for(i in 1:length(glm.obj)){
  nm <- names(glm.obj)[i]
  
  tmp.df <- data.frame("species"=species,
                       "x"=c(1,2) + (i-1.5)/4,
                       "lwr"=rep(0, length(species)),
                       "upr"=rep(0, length(species)),
                       "method"=nm)
  rownames(tmp.df) <- species
  for(s in species){
    tmp.df[s,"y"] <- CalcRateSum(glm.obj[[nm]][[s]], avg.gen[[s]], denom=2*auto.size[s])
    
    tmp <- CalcRateSum(glm.sim[[nm]][[s]], avg.gen[[s]], denom=2*auto.size[s])
    tmp.df[s, c("lwr","upr")] <- quantile(tmp, probs=probs) 
  }
  
  plt.df <- rbind(plt.df, tmp.df)
}

use.xlim <- c(0.5, length(species)+0.5)

pwr <- -8 # Set axes magnitude units
for(j in c('y','lwr','upr')){
  plt.df[[j]] <- plt.df[[j]] / (10^pwr)
}
gg.obj <- ggplot(data=plt.df) + geom_pointrange(aes(x=x, y=y, ymin=lwr, ymax=upr, color=species, shape=method), size=1, fatten=2, fill='white')

gg.obj <- gg.obj + scale_shape_manual(name="Compartment", labels=c("ortho"="Orthologous", "all"="All callable"), values=c("ortho"=21, "all"=22))
gg.obj <- gg.obj + scale_x_continuous(name="", breaks=seq_along(species), labels=species.nms[species], limits=use.xlim)
gg.obj <- gg.obj + scale_color_manual(guide=FALSE, values=sp.color.vec[levels(plt.df$species)])
gg.obj <- gg.obj + ylab(bquote(Mutations ~ per ~ bp%*%10^.(pwr)))
gg.obj <- gg.obj + ggtitle("Mutation rate in and out-\nside of ortho regions")
gg.obj <- gg.obj + bw.theme + axis.theme + xgrid.off.theme
gg.obj <- gg.obj + theme(legend.background=element_rect(fill="white", size=0.5, linetype="solid", color="grey20"),
                         legend.position = c(0.72, 0.84))
print(gg.obj)

out.fn <- JoinPath(res.dir, "figS5A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=3.2, height=4, units="in", useDingbats=FALSE)
}
