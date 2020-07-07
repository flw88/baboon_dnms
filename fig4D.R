#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure 4D

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("ape")
library("adephylo")
library("ggrepel")

source("setup_env.R")
set.seed(1212)


##### LOAD DATA #####
source("load_data.R")

# Number of mutational opportunities of different types 
# (Strong/Weak)
ortho.opp <- fread(JoinPath(dat.dir, "mut.strong_weak.txt"), sep="\t", header=TRUE)
setkey(ortho.opp, TYPE)

##### SETTINGS #####
save.plots <- TRUE


##### FUNCTIONS #####
### Convert a two column data.table to a named vector using nm.col
TwoCol2Vector <- function(tab, nm.col){
  if(ncol(tab) != 2){ stop("tab does not have two columns") }
  val.col <- names(tab)[names(tab) != nm.col]
  out.vec <- tab[, get(val.col)]
  names(out.vec) <- tab[, get(nm.col)]
  return(out.vec)
}

###
CountTypes <- function(x, types=c("CpG S>W", "S>S", "W>S", "W>W", "non-CpG S>W"), 
                       return.df=FALSE){
  z <- table(x)
  if(!all(names(z) %in% types)){ stop("Some elements of x are not of specified types") }

  z[ types[!(types %in% names(z))] ] <- 0
  z <- z[types]

  if(return.df){
    z <- data.frame("type"=names(z), "n"=as.vector(z))
  }

  return(z)
}

CalcRates <- function(count.tab, fdr1, fdr2, fnr, denom){
  TypeToOpp <- function(nm){
  # Convert mutation type to mutational opportunity 
  # ("S>W" becomes "S", "CpG S>W" becomes "CpG")
    return( gsub(nm, pattern="[ |>].+", replacement="") )
  }

  count.tab[CLUST_SIZE == 1, y := n*(1-fdr1[F1]), by=.(F1)]
  count.tab[CLUST_SIZE == 2, y := n*(1-fdr2),     by=.(F1)]
  # out.tab <- count.tab[, .( rate = 2*sum(y)/((1-fnr[F1])*denom[TypeToOpp(type), get(F1)[1]]) ), by=.(F1, type)]
  out.tab <- count.tab[, .( rate = 2*sum(y)/((1-fnr[F1])*denom[TypeToOpp(type), SIZE]) ), by=.(F1, type)]
  out.tab <- out.tab[, .(mu = sum(rate)), by=.(type)]
  setkey(out.tab, type)
  return(out.tab)
}

ResampleBlock <- function(block){
  block.uniq <- unique(block)
  inds <- sample(block.uniq, size=length(block.uniq), replace=TRUE)
  inds <- unname(unlist(sapply(inds, function(x) which(block == x))))
  return(inds) # Returns indices of block array that correspond to the resampled elements
}

### Convert Strong/Weak type to Base types
# (e.g., S>S becomes C>G)
ConvertToSW <- function(base.type, is.cpg){
  if       (base.type == "C>G"){
    sw.type <- "S>S"
  } else if(base.type == "T>A"){
    sw.type <- "W>W"
  } else if(base.type %in% c("T>C", "T>G")){
    sw.type <- "W>S"
  } else if((base.type == "CpG>TpG") || ((base.type == "C>A") &&  is.cpg)){
    sw.type <- "CpG S>W"
  } else if((base.type == "CpH>TpH") || ((base.type == "C>A") && !is.cpg)){
    sw.type <- "non-CpG S>W"
  } else {
    stop("Couldn't convert to SW type")
  }
  return(sw.type)
}


##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]

  # Remove non-orthologous DNMs
  dnm.df[[s]] <- dnm.df[[s]][ORTHO == 1]

  # As is, the TYPE column only distinguishes CpG from non-CpG
  # in the case of C>T, but we also need to check C>A if considering
  # all S>W types
  dnm.df[[s]][, TYPE := ConvertToSW(TYPE, CPG), by=.(CHROM, POS, F1)]
}

##### COMPUTE FDR #####
fdr.ones <- list()
fdr.twos <- list()
for(s in species){
  fdr.ones[[s]] <- ComputeSingletFDR(tri.df[s], dnm.df[[s]])
  fdr.twos[[s]] <- ComputeDoubletFDR(tri.df[s], dnm.df[[s]])
}

# Estimate FDR for baboon F1s without an F2
s <- "bab"
RegressFDR(tri.df[s], fdr.ones[[s]])

##### DOWNLOAD SUBSTIUTITON DATA #####
# Download substitution data from Moorjani et al. GitHub if needed
sub.url <- "https://raw.githubusercontent.com/priyamoorjani/Molecular-clock_figures-and-data/master/data/Figure3b.txt"
sub.fn <- JoinPath(dat.dir, basename(sub.url))
if(!file.exists(sub.fn)){
  download.file(sub.url, sub.fn)
}

# Mutation/substitution types to compute relative rates for
sw.data <- list()
for(s in species){
  sw.data[[s]] <- data.table("type" =c("S>S", "W>W", "W>S", "CpG S>W", "non-CpG S>W"),
                             "type2"=c("ss", "ww", "ws", "cpg_sw", "noncpg_sw"), # Names used by Moorjani et al.
                             "k"    =as.double(NA), # Substution rates per yr
                             "mu"   =as.double(NA)) # Mutation rates per yr
  setkey(sw.data[[s]], type)
}

# Obtain substitution rates and trees
for(nm in sw.data[["hum"]][,type]){
  tmp <- ComputeSubRates(sub.fn, sw.data[["hum"]][nm, type2])
  sw.data[["hum"]][nm, k := tmp[["rates"]]["hum"] ]
  sw.data[["bab"]][nm, k := tmp[["rates"]]["bab"] ]
}


##### CALCULATE MUTATION RATES #####
# (Also start building dataframe for plotting)
plt.df <- data.table("type"    =sw.data[["hum"]][,type],
                     "k.ratio" =as.double(NA),
                     "mu.ratio"=as.double(NA),
                     "mu.ratio.lwr"=as.double(NA),
                     "mu.ratio.upr"=as.double(NA))
setkey(plt.df, type)

n.sim <- 1000
u.sim <- list()
type.nms <- plt.df[,type]
for(s in species){
  rm.f1s <- tri.df[s][is.na(Mat_Age), unique(F1)]
  t.df <- tri.df[s][!(F1 %in% rm.f1s)]
  # Remove F1s with no recorded Maternal Age at Reproduction
  d.df <- dnm.df[[s]][!(F1 %in% rm.f1s)]

  count.tab <- d.df[, CountTypes(TYPE, return.df=TRUE), by=.(F1, CLUST_SIZE)]

  fdr1  <- TwoCol2Vector(fdr.ones[[s]][, .(F1, est)], nm.col="F1")
  fdr2  <- fdr.twos[[s]]
  fnr   <- TwoCol2Vector(t.df[!duplicated(F1), .(F1,FNR)], nm.col="F1")
  # denom <- ortho.opp[, c("TYPE", t.df[!(duplicated(F1)), F1]), with=FALSE]
  denom <- ortho.opp[, .(TYPE, SIZE = get(s) * 2)]
  

  age.sum <- t.df[ !duplicated(F1), sum(Pat_Age+Mat_Age-(2*avg.con[s])) ]

  u <- CalcRates(count.tab, fdr1, fdr2, fnr, denom)
  sw.data[[s]][, mu := u[type, mu / age.sum]]

  # Resample by block for CIs
  u.sim[[s]] <- matrix(0.0, nrow=length(type.nms), ncol=n.sim)
  rownames(u.sim[[s]]) <- type.nms
  
  for(i in 1:n.sim){
    r <- d.df[, .I[ResampleBlock(BLOCK)], by=F1]$V1 # Row indices
    count.tab <- d.df[r, CountTypes(TYPE, return.df=TRUE), by=.(F1, CLUST_SIZE)]
    u <- CalcRates(count.tab, fdr1, fdr2, fnr, denom)
    u.sim[[s]][type.nms,i] <- u[type.nms, mu / age.sum]
  }

}

# Compute point estimate of ratio
plt.df[type.nms, k.ratio := sw.data[["bab"]][type.nms,k] / sw.data[["hum"]][type.nms,k] ]
plt.df[type.nms, mu.ratio:= sw.data[["bab"]][type.nms,mu]/ sw.data[["hum"]][type.nms,mu]]

# Confidence intervals
q <- t( apply(u.sim[["bab"]] /u.sim[["hum"]], 1, function(x) quantile(x, probs=c(0.025, 0.975))) )
plt.df[type.nms, mu.ratio.lwr := q[type.nms,1]]
plt.df[type.nms, mu.ratio.upr := q[type.nms,2]]


##### PLOT RELATIVE RATES SUB VS MUT #####
# Add columns for BGC status, color
plt.df[c("S>S", "W>W"),    c("BGC.status", "color") := list("not sensitive", "#1C7FDB")] # Blue
plt.df[grepl("S>W", type), c("BGC.status", "color") := list("sensitive (-)", "#C93203")] # Dark red
plt.df["W>S",              c("BGC.status", "color") := list("sensitive (+)", "#FF7144")] # Light red


plt.df[, nudge.x := 0.0]
plt.df["S>S", nudge.x := -0.15]
plt.df["W>W", nudge.x := -0.10]

plt.df[, nudge.y := 0.0]
plt.df[        "S>S", nudge.y := -0.10]
plt.df[        "W>W", nudge.y :=  0.20]
# plt.df["non-CpG S>W", nudge.y := -0.25]
plt.df["CpG S>W", nudge.y := -0.10]

use.colors <- TwoCol2Vector(plt.df[,.(BGC.status, color)], nm.col="BGC.status")
plt.lim <- c(0.5,2) # Plot limits

gg.obj <- ggplot(aes(x=mu.ratio, y=k.ratio, label=type, color=BGC.status, 
                     xmin=mu.ratio.lwr, xmax=mu.ratio.upr), data=plt.df)
gg.obj <- gg.obj + geom_abline(slope=1, intercept=0, alpha=0.2)
gg.obj <- gg.obj + geom_errorbarh(height=0, size=0.75) 
gg.obj <- gg.obj + geom_point(fill="white", shape=21, size=2, stroke=0.75)
gg.obj <- gg.obj + geom_text_repel(show.legend=FALSE, size=3, 
                                   nudge_x=plt.df$nudge.x, nudge_y=plt.df$nudge.y, seed=2)
gg.obj <- gg.obj + scale_color_manual(values=use.colors, name="BGC status") # Use parent colors
gg.obj <- gg.obj + xlab("Relative Yearly Mutation Rate") + ylab("Relative Yearly Substitution Rate")
gg.obj <- gg.obj + ggtitle("Comparison of Relative (Baboon/\nHuman) Substitution and Mutation\nRates by Type")
gg.obj <- gg.obj + coord_cartesian(xlim=plt.lim, ylim=plt.lim)
gg.obj <- gg.obj + bw.theme + axis.theme
gg.obj <- gg.obj + theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="grey20"),
                         legend.position = c(0.2,0.19))
print(gg.obj)

out.fn <- JoinPath(res.dir, "fig4D.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=4, height=4.4, units="in", useDingbats=FALSE)
}
