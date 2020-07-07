#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure SB

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("arm")
library("digest")

source("setup_env.R")
set.seed(8080)


##### LOAD DATA #####
source("load_data.R")


##### SETTINGS #####
save.plots <- TRUE


##### FUNCTIONS #####
### Builds the dataframe object from which the regression
### plot is made
BuildRegrPlotDf <- function(cur.glm, ci=0.95, length.out=100){
  range.age <- range(cur.glm$data$age)
  out.df <- data.table("age" = SeqRange(range.age, length.out=length.out))
  tmp <- predict(cur.glm, newdata=out.df, se.fit=T)
  out.df[, dnm := tmp$fit ]
  out.df[, lwr := dnm + qnorm((1-ci)/2) * tmp$se.fit ]
  out.df[, upr := dnm + qnorm((1+ci)/2) * tmp$se.fit ]
  return(out.df)
}

##### LOAD REGRESSION RESULTS #####
glm.fn <- JoinPath(res.dir, "glm.ortho.RData")
StopNoGLMData(glm.fn)

load(glm.fn) # Loads save.list
glm.obj <- save.list$glm.obj
glm.sim <- save.list$glm.sim


### Using age post-puberty
glm.pub.fn <- JoinPath(res.dir, "glm.ortho.postpub.RData")
StopNoGLMData(glm.pub.fn)

load(glm.pub.fn) # Loads save.list
glm.pub <- save.list$glm.obj
glm.pub.sim <- save.list$glm.sim


##### PLOT SPECIES TOGETHER #####
pts.df <- NULL; lines.df <- NULL
for(s in species){
  for(p in parents){
    tmp.df <- data.table::copy(glm.obj[[s]][[p]]$data[, .(age, dnm)])
    tmp.df[, c("Species", "Parent") := .(species.nms[s], paste0(p, "ernal"))]
    pts.df <- rbindlist(list(pts.df , tmp.df))
    
    tmp.df <- BuildRegrPlotDf(glm.obj[[s]][[p]])
    tmp.df[, c("Species", "Parent") := .(species.nms[s], paste0(p, "ernal"))]
    lines.df <- rbindlist(list(lines.df, tmp.df))
  }
}

pts.df[,   GRP := paste(Species, Parent, sep="_")]
lines.df[, GRP := paste(Species, Parent, sep="_")]

shp.vals <- c("Paternal"=21, "Maternal"=24)
fill.vals <- c("Human"=unname(sp.color.vec["hum"]), 
               "Baboon"=unname(sp.color.vec["bab"]))
col.vals <- fill.vals
lin.vals <- c("Paternal"=1, "Maternal"=2)



gg.obj <- ggplot()
gg.obj <- gg.obj + geom_line(aes(x=age, y=dnm, color=Species, linetype=Parent, group=GRP), data=lines.df, alpha=1, size=0.75)
gg.obj <- gg.obj + geom_ribbon(aes(x=age, ymin=lwr, ymax=upr, fill=Species, group=GRP), data=lines.df, alpha=0.3)
gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm, fill=Species, shape=Parent), data=pts.df, color="#333333", size=1.8, stroke=0.8) 
gg.obj <- gg.obj + scale_shape_manual(values=shp.vals, guide=FALSE) + scale_linetype_manual(values=lin.vals, guide=FALSE)
gg.obj <- gg.obj + scale_fill_manual(values=fill.vals, guide=FALSE) + scale_color_manual(   values=col.vals, guide=FALSE)
gg.obj <- gg.obj + bw.theme + axis.theme
gg.obj <- gg.obj + xlab("Parental Age at Conception (y)") + ylab("Inferred Number of DNMs")
print(gg.obj)

# Save plot
out.fn <- JoinPath(res.dir, "figS7A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=8, height=4, units="in", useDingbats=FALSE)
}

##### REGRESS EACH SPECIES USING AGE PUBERTY (NOT CONCEPTION) #####

### Plot
pts.df <- NULL; lines.df <- NULL
p <- "Pat"
for(s in species){
  tmp.df <- data.table::copy(glm.pub[[s]][[p]]$data[, .(age, dnm)])
  tmp.df[, c("Species", "Parent") := .(species.nms[s], paste0(p, "ernal"))]
  pts.df <- rbindlist(list(pts.df , tmp.df))
  
  tmp.df <- BuildRegrPlotDf(glm.pub[[s]][[p]])
  tmp.df[, c("Species", "Parent") := .(species.nms[s], paste0(p, "ernal"))]
  lines.df <- rbindlist(list(lines.df, tmp.df))
}

pts.df[,   GRP := paste(Species, Parent, sep="_")]
lines.df[, GRP := paste(Species, Parent, sep="_")]

shp.vals <- c("Paternal"=21, "Maternal"=24)
fill.vals <- c("Human"=unname(sp.color.vec["hum"]), 
               "Baboon"=unname(sp.color.vec["bab"]))

col.vals <- fill.vals
lin.vals <- c("Paternal"=1, "Maternal"=2)


gg.obj <- ggplot()
gg.obj <- gg.obj + geom_line(aes(x=age, y=dnm, color=Species, linetype=Parent, group=GRP), data=lines.df, alpha=1, size=0.75)
gg.obj <- gg.obj + geom_ribbon(aes(x=age, ymin=lwr, ymax=upr, fill=Species, group=GRP), data=lines.df, alpha=0.3)
gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm, fill=Species, shape=Parent), data=pts.df, color="#333333", size=1.8, stroke=0.8) 
gg.obj <- gg.obj + scale_shape_manual(values=shp.vals, guide=FALSE) + scale_linetype_manual(values=lin.vals, guide=FALSE)
gg.obj <- gg.obj + scale_fill_manual(values=fill.vals, guide=FALSE) + scale_color_manual(   values=col.vals, guide=FALSE)
gg.obj <- gg.obj + bw.theme + axis.theme
gg.obj <- gg.obj + xlab("Paternal Years Post-Puberty (at Conception)") + ylab("Inferred Number of DNMs")
print(gg.obj)

# Save plot
out.fn <- JoinPath(res.dir, "figS7B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=6, height=4, units="in", useDingbats=FALSE)
}

##### VARIANCE/DEVIANCE EXPLAINED BY PATERNAL AGE #####
VarianceExplained <- function(cur.lm){
  return(summary(cur.lm)$r.squared)
}

DevianceExplained <- function(cur.glm){
  return(1 - cur.glm$deviance / cur.glm$null.deviance)
}

cat("###################################################\n",
    "########### Deviance/Variance Explained ###########\n", sep="")
lm.df <- list()
for(s in species){
  lm.df[[s]] <- list()
  for(p in parents){
    lm.df[[s]][[p]] <- data.table::copy(glm.obj[[s]][[p]]$data)
  }
}

lm.obj <- list()
for(s in species){
  lm.obj[[s]] <- list()
  for(p in parents){
    lm.obj[[s]][[p]] <- lm(dnm ~ age, data=lm.df[[s]][[p]])
    cat(sprintf("=====%7s , %sernal =====\n", species.nms[s], p))
    
    var.expl <- VarianceExplained(lm.obj[[s]][[p]])
    cat(sprintf("%20s: %2.1f%%\n", "Variance explained", 100*var.expl))
    
    dev.expl <- DevianceExplained(glm.obj[[s]][[p]])
    cat(sprintf("%20s: %2.1f%%\n", "Deviance explained", 100*dev.expl))
  }
}

