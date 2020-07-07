#!/usr/bin/env Rscript
rm(list=ls())

# Script for

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("digest")

source("setup_env.R")
set.seed(1212)


##### LOAD DATA #####
source("load_data.R")

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
TabulateDNMs <- function(x, use.levels=unique(x), return.df=FALSE){
  z <- table(factor(x, levels = use.levels))

  if(return.df){
    z <- data.frame("bin"=names(z), "n"=as.vector(z))
  }

  return(z)
}

CalcRates <- function(count.tab, fdr1, fdr2, fnr, denom){
  count.tab[CLUST_SIZE == 1, y := n*(1-fdr1[F1]), by=.(F1)]
  count.tab[CLUST_SIZE == 2, y := n*(1-fdr2),     by=.(F1)]
  
  out.tab <- count.tab[, .( rate = sum(y)/((1-fnr[F1])*as.numeric(denom[bin, get(F1)[1]])) ), by=.(F1, bin)]
  out.tab <- out.tab[, .(mu = mean(rate)), by=.(bin)]
  setkey(out.tab, bin)
  return(out.tab)
}

BinPairName <- function(bin.pair){
  return( paste(bin.pair, collapse="v") )
}

CalcPval <- function(count.tab, denom, bin.pair.mat=rbind(c("0", "1"))){
  # Conditional test for Poisson rates Xi ~ Pois
  # X1 + X2 = k then X1 | k ~ Binomial
  # NB: assumes that underlying Poisson rate parameter is the same for each individual
  # (i.e., ignores age effect)
  # bin.pair.mat indicates which bins to compare
  
  p.val.tab <- data.table("Comparison" = apply(bin.pair.mat, 1, BinPairName),
                          "PVAL" = as.numeric(NA))
  setkey(p.val.tab, Comparison)
  for(i in 1:nrow(bin.pair.mat)){
    bin.pair <- bin.pair.mat[i,]
    k <- c(); n <- c()
    for(b in bin.pair){
      k[b] <- count.tab[bin == b, sum(n)]
      n[b] <- rowSums(denom[BIN == b, -"BIN"])
    }
    k.sum <- sum(k)
    
    # Probability of binomial distribution under assumption that rates are the same
    binom.prob <- (n[1]/n[2]) / (1 + n[1]/n[2])
    
    p <- sum(dbinom(x=(k[1]:k.sum), size=k.sum, prob=binom.prob))
    p.val.tab[BinPairName(bin.pair), PVAL := p]
  }
  
  return(p.val.tab)
}

ResampleBlock <- function(block){
  block.uniq <- unique(block)
  inds <- sample(block.uniq, size=length(block.uniq), replace=TRUE)
  inds <- unname(unlist(sapply(inds, function(x) which(block == x))))
  return(inds) # Returns indices of block array that correspond to the resampled elements
}


##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]
  dnm.df[[s]][, "ORTHO.REPEAT" := paste(ORTHO, REPEAT, sep="")]
}

##### GET ORTHOLOGOUS AND NON-ORTHOLOGOUS CALLABLE GENOME SIZES #####
clb.gnm <- list()
for(s in species){
  clb.gnm[[s]] <- list()
  clb.gnm[[s]][["ORTHO"]]  <- data.table(BIN=c("0", "1"))
  clb.gnm[[s]][["REPEAT"]] <- data.table(BIN=c("0", "1"))
  clb.gnm[[s]][["ORTHO.REPEAT"]] <- data.table(BIN=c("00", "01", "10", "11"))
  for(f1 in unique(tri.df[s, F1])){
    cur.clb <- tri.df[F1 == f1, Clb][1]
    cur.clb.rep <- tri.df[F1 == f1, Clb_Rep][1]
    cur.ortho.rep <- unname(rep.size$ortho[s])
    
    new.val <- 2 * c(cur.clb - ortho.size, ortho.size)
    clb.gnm[[s]][["ORTHO"]][, eval(f1) := new.val]
    
    new.val <- 2 * c(cur.clb - cur.clb.rep, cur.clb.rep)
    clb.gnm[[s]][["REPEAT"]][, eval(f1) := new.val]
    
    new.val <- 2 * c(cur.clb - cur.clb.rep - ortho.size + cur.ortho.rep, 
                     cur.clb.rep - cur.ortho.rep, 
                     ortho.size - cur.ortho.rep, 
                     rep.size$ortho[s])
    clb.gnm[[s]][["ORTHO.REPEAT"]][, eval(f1) := new.val]
  }
  
  # Set keys
  for(cmp in names(clb.gnm[[s]])){
    setkey(clb.gnm[[s]][[cmp]], BIN)
  }
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


##### CALCULATE MUTATION RATES, P-VALUES #####
compartments <- names(clb.gnm$hum)
n.sim <- 1000

test.comparisons <- list("ORTHO"  = rbind(c("0", "1")),
                         "REPEAT" = rbind(c("1", "0")),
                         "ORTHO.REPEAT" = rbind(c("00", "10"),
                                                c("00", "01"),
                                                c("10", "11")))

# Memoize since it takes a while to calculate u.sim
dig.str <- digest(list(tri.df, dnm.df, clb.gnm, n.sim), algo="md5")
mem.fn <- JoinPath(res.dir, sprintf("repeat2.%s.RData", dig.str))
if(file.exists(mem.fn)){
  cat("Loading u and u.sim from", mem.fn, "\n")
  load(mem.fn)
  u     <- mem.list[["u"]]
  u.sim <- mem.list[["u.sim"]]
  p.val <- mem.list[["p.val"]]
} else {
  u <- list()
  u.sim <- list()
  p.val <- list()
  
  for(s in species){
    u[[s]] <- list()
    u.sim[[s]] <- list()
    p.val[[s]] <- list()
    for(cmp in compartments){ # Iterate over comparisons
      t.df <- tri.df[s]
      d.df <- dnm.df[[s]]
      
      denom <- clb.gnm[[s]][[cmp]]
      
      count.tab <- d.df[, TabulateDNMs(as.character(get(cmp)), use.levels=denom[["BIN"]], return.df=TRUE), by=.(F1, CLUST_SIZE)]
      count.tab[, bin := as.character(bin)]
    
      fdr1  <- TwoCol2Vector(fdr.ones[[s]][, .(F1, est)], nm.col="F1")
      fdr2  <- fdr.twos[[s]]
      fnr   <- TwoCol2Vector(t.df[!duplicated(F1), .(F1,FNR)], nm.col="F1")
      
      p.val[[s]][[cmp]] <- CalcPval(count.tab, denom, bin.pair.mat=test.comparisons[[cmp]])
    
      u[[s]][[cmp]] <- CalcRates(count.tab, fdr1, fdr2, fnr, denom)
      
      # Resample by block for CIs
      u.sim[[s]][[cmp]] <- matrix(0.0, nrow=denom[,.N], ncol=n.sim)
      rownames(u.sim[[s]][[cmp]]) <- denom[["BIN"]]
      
      for(i in 1:n.sim){
        r <- d.df[, .I[ResampleBlock(BLOCK)], by=F1]$V1 # Row indices
        count.tab <- d.df[r, TabulateDNMs(as.character(get(cmp)), use.levels=denom[["BIN"]], return.df=TRUE), by=.(F1, CLUST_SIZE)]
        count.tab[, bin := as.character(bin)]
        u.sim[[s]][[cmp]][,i] <- CalcRates(count.tab, fdr1, fdr2, fnr, denom)[, mu]
      }
      
    }
  
  }
  
  mem.list <- list("u"=u, "u.sim"=u.sim, "p.val"=p.val)
  cat("Saving u and u.sim to", mem.fn, "\n")
  save(mem.list, file=mem.fn)
}

##### BUILD RESULTS TABLE #####
res.df <- list() 
for(cmp in compartments){
  res.df[[cmp]] <- NULL
  
  for(s in species){
    q <- t( apply(u.sim[[s]][[cmp]], 1, function(x) quantile(x, probs=c(0.025, 0.975))) )
    cur.df <- data.table("BIN"    = u[[s]][[cmp]][, bin],
                         "MU"     = u[[s]][[cmp]][, mu ],
                         "MU_LWR" = q[,1],
                         "MU_UPR" = q[,2],
                         "Species"= species.nms[s],
                         "Sps"    = s)
    res.df[[cmp]] <- rbindlist(list(res.df[[cmp]], cur.df))
  }
  
  res.df[[cmp]][, Sps := factor(Sps, levels=species)]
  res.df[[cmp]][, Species := factor(Species, levels=unname(species.nms))]
  setkey(res.df[[cmp]], Sps, BIN)
}

##### PLOT RESULTS #####

use.colors <- sp.color.vec[species]
names(use.colors) <- species.nms[species]
use.fill <- c("0"="#FFFFFF", "1"="#444444")
use.shapes <- c("0"=21, "1"=22)


pwr <- -8 # Set axes magnitude units
cmp <- "ORTHO.REPEAT"
plt.df <- data.table::copy(res.df[[cmp]])
plt.df[, ":=" (x  = as.numeric(.I),
               MU = MU / (10^pwr),
               MU_LWR = MU_LWR / (10^pwr),
               MU_UPR = MU_UPR / (10^pwr))]
plt.df[Sps == "bab", x := x + 0.5]
x.vline <- mean(c(plt.df[Sps == "hum", max(x)], plt.df[Sps == "bab", min(x)]))
x.breaks <- c(mean(plt.df[Sps=="hum", range(x)]), mean(plt.df[Sps=="bab", range(x)]))

gg.obj <- ggplot(data=plt.df)

# plt.df[, c("ORTHO", "REPEAT") := tstrsplit(BIN, split="") ]
# gg.obj <- gg.obj + geom_pointrange(aes(x=x, y=MU, ymin=MU_LWR, ymax=MU_UPR, color=Species, shape=ORTHO, fill=REPEAT), size=1, fatten=2)
gg.obj <- gg.obj + geom_pointrange(aes(x=x, y=MU, ymin=MU_LWR, ymax=MU_UPR, color=Species, shape=BIN, fill=BIN), size=1, fatten=2)

cur.labs <- c("Non-orthologous, non-repetitive", "Non-orthologous, repetitive", "Orthologous, non-repetitive", "Orthologous, repetitive")
names(cur.labs)   <- c("00", "01", "10", "11")
cur.shapes <- c(use.shapes["0"], use.shapes["0"], use.shapes["1"], use.shapes["1"])
cur.fill   <- c(use.fill["0"],   use.fill["1"],   use.fill["0"],   use.fill["1"])
names(cur.shapes) <- names(cur.labs)
names(cur.fill)   <- names(cur.labs)

gg.obj <- gg.obj + scale_shape_manual(name="Compartment", labels=cur.labs, values=cur.shapes)
gg.obj <- gg.obj + scale_fill_manual(name="Compartment", labels=cur.labs, values=cur.fill)

gg.obj <- gg.obj + geom_vline(xintercept=x.vline, color="gray")
gg.obj <- gg.obj + scale_color_manual(guide=FALSE, values=use.colors)
gg.obj <- gg.obj + scale_x_continuous(name="", breaks=x.breaks, labels=species.nms, limits=c(-0.5,0.5) + plt.df[, range(x)])
gg.obj <- gg.obj + ylab(bquote(Mutations ~ per ~ bp%*%10^.(pwr)))
gg.obj <- gg.obj + ggtitle(cmp)
gg.obj <- gg.obj + bw.theme + axis.theme + xgrid.off.theme
gg.obj <- gg.obj + theme(legend.background=element_rect(fill="white", size=0.5, linetype="solid", color="grey20"), 
                         legend.position = c(0.68, 0.78))
print(gg.obj)

out.fn <- JoinPath(res.dir, "figS5B.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=4.8, height=4, units="in", useDingbats=FALSE)
}

##### PRINT REPORT #####
if(save.plots){
  out.fn <- JoinPath(res.dir, "report.ortho_repeat.mut_rate.tsv")
} else {
  out.fn <- ""
}

cat("", append=FALSE, file=out.fn)
for(s in species){
  cat("###################################\n", append=TRUE, file=out.fn)
  cat(species.nms[s], "\n", append=TRUE, file=out.fn)
  cat("##########\n", append=TRUE, file=out.fn)
  for(cmp in compartments){
    cat(sprintf("===== %s =====\n", cmp), append=TRUE, file=out.fn)
    for(i in 1:nrow(test.comparisons[[cmp]])){
      bin.pair <- test.comparisons[[cmp]][i,]
      
      ratio <- u[[s]][[cmp]][bin.pair[1], mu] / u[[s]][[cmp]][bin.pair[2], mu]
      ratio.sim <- u.sim[[s]][[cmp]][bin.pair[1],] / u.sim[[s]][[cmp]][bin.pair[2],]
      
      nm <- BinPairName(bin.pair)
      p <-p.val[[s]][[cmp]][nm, PVAL]
      
      if(ratio < 1){
        ratio <- 1 / ratio
        ratio.sim <- 1 / ratio.sim
        nm <- BinPairName(rev(bin.pair))
        p <- 1 - p
      }
      
      ci <- quantile(ratio.sim, probs=c(0.025, 0.975))
      tmp.str <- sprintf("  %s :   ratio = %0.2f (%0.2f -- %0.2f) , p-value = %0.2e \n", nm, ratio, ci[1], ci[2], p)
      cat(tmp.str, append=TRUE, file=out.fn)
    }
    
    cat("\n", append=TRUE, file=out.fn)
  }
}

