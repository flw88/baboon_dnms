#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure 4

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("digest")
library("ape")
library("adephylo")
library("ggtree")
library("ggrepel")

source("setup_env.R")


##### LOAD DATA #####
source("load_data.R")

##### SETTINGS #####
save.plots <- TRUE


##### LIFE HISTORY TEST VALUES #####

# Range of paternal generation times
pat.gen.range <- list('hum'=c(3, avg.gen$hum["Pat"]),
                      'bab'=c(3, 15))
for(s in species){ names(pat.gen.range[[s]]) <- c("min", "max") }

# Range of pat / mat generation time ratios
gen.ratio <- list('hum'=seq(0.8, 1.2, by=0.1),
                  'bab'=seq(0.8, 1.2, by=0.1))

##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]

  # Remove non-orthologous DNMs
  dnm.df[[s]] <- dnm.df[[s]][ORTHO == 1]
}

##### DOWNLOAD SUBSTIUTITON DATA #####
sub.url <- "https://raw.githubusercontent.com/priyamoorjani/Molecular-clock_figures-and-data/master/data/Figure1.txt"
sub.fn <- JoinPath(dat.dir, basename(sub.url))
if(!file.exists(sub.fn)){
  download.file(sub.url, sub.fn)
}

tmp <- ComputeSubRates(sub.fn, "average_rate")
k.y <- tmp[["rates"]]
sub.tree <- tmp[["tree"]]


##### LOAD REGRESSION RESULTS #####
glm.fn <- JoinPath(res.dir, "glm.ortho.RData")
StopNoGLMData(glm.fn)

load(glm.fn) # Loads save.list
glm.obj <- save.list$glm.obj
glm.sim <- save.list$glm.sim


##### PLOT TREE WITH OUTGROUP #####
tree.sp <- c("human","baboon","marmoset")
plt.tree <- ComputeSubRates(sub.fn, "average_rate", kp.species=tree.sp)[["tree"]]


gg.obj <- ggtree(plt.tree, ladderize=FALSE, color="#333333", size=1) + theme_tree2() + geom_tiplab(color=c(sp.color.vec, "black"), offset=.001)
# Add relative baboon/human rate to tree
y.pos <- 2.1
line.df <- data.frame("x"=distRoot(plt.tree, "baboon") - Dist2MRCA(plt.tree, c("human","baboon"))["baboon"], 
                      "xend"=distRoot(plt.tree, "baboon"),
                      "y"=y.pos, "yend"=y.pos)
line.col <- "#5738B4"
gg.obj <- gg.obj + geom_segment(aes(x=x, xend=xend, y=y, yend=yend), data=line.df, color=line.col,
                                arrow=arrow(length=unit(0.03, "npc"), ends="both", type="closed"))
gg.obj <- gg.obj + geom_text(aes(x=mean(c(line.df$x, line.df$xend)), y=y.pos, 
                                 label=sprintf("%0.2f", k.y["bab"]/k.y["hum"])), 
                             nudge_y=0.1, hjust="center", color=line.col)

gg.obj <- gg.obj + xlab("Substitutions per bp") + xlim(0, max(plt.tree$edge.length)+0.02)
gg.obj <- gg.obj + theme(axis.line.x=element_line(), axis.text=element_text(size=11))
gg.obj <- gg.obj + ggtitle("All Mutation Types")
print(gg.obj)

out.fn <- JoinPath(res.dir, "fig4A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=5, height=2.9, units="in", useDingbats=FALSE)
}

##### SEX AVERAGED RATE PER YEAR #####
CalcRatePerYear <- function(obj.list, age.list, denom){ 
  # Takes in list of "glm" or "sim" objects and ages (names=c("Pat", "Mat"))
  mu.g <- 0
  for(p in names(obj.list)){
    mu.g <- mu.g + CalcRateAtAge(coef(obj.list[[p]]), age.list[[p]], denom)
  }
    
  mu.y <- 2 * mu.g / sum(age.list)
  return(mu.y) # Return mutation rates per year for sampled coefficients
}

ci <- 0.95; probs <- c(1-ci, 1+ci) / 2
mu.g <- list() # Point estimate of mutation rate per year
mu.y <- c(); mu.y.sim <- list(); mu.y.ci <- list()
for(s in species){
  mu.g[[s]] <- c()
  for(p in parents){
    mu.g[[s]][p] <- CalcRateAtAge(glm.obj[[s]][[p]], avg.gen[[s]][p], denom=2*auto.size[s])
  }
  mu.y[s] <- 2*sum(mu.g[[s]]) / sum(avg.gen[[s]])
  
  mu.y.sim[[s]] <- CalcRatePerYear(glm.sim[[s]], avg.gen[[s]], denom=2*auto.size[s])
  
  mu.y.ci[[s]] <- quantile(mu.y.sim[[s]], probs=probs)
  names(mu.y.ci[[s]]) <- c("lwr", "upr")
}


##### RELATIVE RATE PER YR #####
rel.y <- mu.y["bab"] / mu.y["hum"]
rel.y.ci <- quantile(mu.y.sim$bab / mu.y.sim$hum, probs=probs)
names(rel.y.ci) <- c("lwr","upr")

# Report
cat(sprintf("Relative rates [bab/hum]: %0.2f (%d%% CI: %0.2f, %0.2f)\n",
            mean(rel.y), ci*100, rel.y.ci["lwr"], rel.y.ci["upr"]))


##### PLOT PER-YEAR RATE #####
plt.df <- data.table("species"=species,
                     "x"  =seq_along(species),
                     "y"  =rep(0.0, length(species)),
                     "lwr"=rep(0.0, length(species)),
                     "upr"=rep(0.0, length(species)))
setkey(plt.df, species)
use.xlim <- c(0.5, nrow(plt.df)+0.5)

for(s in species){
  plt.df[s, y := mu.y[s]]

  plt.df[s, lwr := mu.y.ci[[s]]["lwr"]]
  plt.df[s, upr := mu.y.ci[[s]]["upr"]]
}

pwr <- -9 # Set axes magnitude units
for(j in c("y","lwr","upr")){
  plt.df[[j]] <- plt.df[[j]] / (10^pwr)
}

gg.obj <- ggplot(data=plt.df) + geom_pointrange(aes(x=x, y=y, ymin=lwr, ymax=upr, color=species), shape=21, size=1, fatten=2, fill="white")
gg.obj <- gg.obj + scale_x_continuous(name="", breaks=plt.df$x, labels=species.nms[plt.df$species], limits=use.xlim)
gg.obj <- gg.obj + scale_color_manual(values=sp.color.vec[plt.df$species])
gg.obj <- gg.obj + ylab(bquote(Mutations ~ per ~ bp%*%10^.(pwr)))
gg.obj <- gg.obj + ggtitle(sprintf("Mutation rate per year\n(using mean population age)"))
gg.obj <- gg.obj + bw.theme + axis.theme + xgrid.off.theme
print(gg.obj)

out.fn <- JoinPath(res.dir, "fig4B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=4, height=4, units="in", useDingbats=FALSE)
}


##### DIVERGENCE TIME USING DECODE AGE EFFECTS #####
fossil.bound <- 34
# Use either regression estimates from this study (fig 4C) or from Jonsson et al./Gao et al. (fig S8)
for(bab.mut.param in c("est", "dec")){
  # Vary paternal generation time and gen time ratio
  div.df <- data.frame()
  n.pts <- 100
  for(s in species){
    age.extrema <- range(gen.ratio[[s]])
    names(age.extrema) <- c("min","max")
  
    pat.age <- SeqRange(pat.gen.range[[s]], length.out=n.pts)
    tmp.df <- data.frame( "x"=pat.age, "species"=as.character(tolower(species.nms[s])) )
    for(j in names(age.extrema)){
      tmp.df[[j]] <- 0
      rat <- age.extrema[j]
      
      use.age <- rbind(pat.age, pat.age/rat)
      rownames(use.age) <- parents
      
      if(s == "hum" || bab.mut.param == "dec"){
        cur.mu.y <- 2*( CalcRateAtAge(gao.est[,"Mat"], age=use.age["Mat",], denom=2*jon.size) + 
                        CalcRateAtAge(gao.est[,"Pat"], age=use.age["Pat",], denom=2*jon.size) )
      } else if(s == "bab" && bab.mut.param == "est"){
        cur.mu.y <- 2*( CalcRateAtAge(coef(glm.obj[[s]]$Mat), age=use.age["Mat",], denom=2*auto.size[s]) + 
                        CalcRateAtAge(coef(glm.obj[[s]]$Pat), age=use.age["Pat",], denom=2*auto.size[s]) )
      }
      cur.mu.y <- cur.mu.y / colSums(use.age)
      tmp.df[[j]] <- (k.y[s] / cur.mu.y) / 1e6
    }
    div.df <- rbind(div.df, tmp.df)
  }
  
  
  # Plot
  color.ordered <- sp.color.vec[substr(levels(div.df$species), 1,3)]
  names(color.ordered) <- NULL
  gg.obj <- ggplot() + geom_ribbon(aes(x=x, ymin=min, ymax=max, fill=species), data=div.df, alpha=0.5)
  gg.obj <- gg.obj + geom_line(aes(x=x, y=min, color=species), data=div.df, size=1) + geom_line(aes(x=x, y=max, color=species), data=div.df, size=1)
  gg.obj <- gg.obj + geom_hline(yintercept = fossil.bound, color="#676767", linetype=2)
  gg.obj <- gg.obj + scale_color_manual(values= c(color.ordered,"grey28"))
  gg.obj <- gg.obj + scale_fill_manual(values=color.ordered)
  cur.row <- div.df[2*n.pts,]
  for(j in names(age.extrema)){
    gg.obj <- gg.obj + geom_text( aes_string(x="x", y=j), data=cur.row, 
                                  label=sprintf("%0.1f", age.extrema[j]),
                                  nudge_x=0.8, color="#5738B4", fontface="bold")
  }
  
  use.title <- "Varying life history traits using\n" 
  if(bab.mut.param == "dec"){
    use.title <- paste0(use.title, "deCODE parameters for both")
  } else if(bab.mut.param == "est") {
    use.title <- paste0(use.title, "deCODE (hum) and est. (bab) parameters")
  }
  
  gg.obj <- gg.obj + ylab("Divergence Time (My)") + xlab("Paternal Age (y)")
  gg.obj <- gg.obj + ggtitle(use.title) + bw.theme + axis.theme
  gg.obj <- gg.obj + theme(legend.background=element_rect(fill="white", size=0.5, linetype="solid", color="grey20"),
                           legend.position = c(0.88, 0.2))
  print(gg.obj)
  
  if(bab.mut.param == "dec"){
    out.fn <- JoinPath(res.dir, "fig4C.pdf")
  } else {
    out.fn <- JoinPath(res.dir, "figS8.pdf")
  }
  if(save.plots){
    cat("Saving", out.fn, "\n")
    ggsave(gg.obj, filename=out.fn, width=5, height=3.5, units="in", useDingbats=FALSE)
  }
}


