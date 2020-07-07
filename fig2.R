#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure 2

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("arm")
library("bbmle")
library("digest")

source("setup_env.R")

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

### Plot sex-specific age effects using ggplot
GGplotAgeEffects <- function(glm.list, use.title="", color.f1s=TRUE, 
                             ref.coef=NULL, fill.col=NULL, stroke.col=NULL){
  # Default fill colors
  if(is.null(fill.col)){
    fill.col <- c("#0000B2", "#B20000")
    names(fill.col) <- parents
  } else if(!all(names(fill.col) == parents) || (length(fill.col) != 2)){
    stop("Bad fill.col")  
  }
  
  # Default stroke colors
  if(is.null(stroke.col)){
    stroke.col <- c("#00007F", "#7F0000")
    names(stroke.col) <- parents
  } else if(!all(names(stroke.col) == parents) || (length(stroke.col) != 2)){
    stop("Bad stroke.col")
  }
  
  
  # Round if need be
  pts.df <- list()
  for(p in names(glm.list)){
    pts.df[[p]] <- as.data.table(glm.list[[p]]$data)
  }
  
  # Build table for plotting fit lines
  lines.df <- list()
  for(p in names(glm.list)){
    lines.df[[p]] <- BuildRegrPlotDf(glm.list[[p]])
  }
  
  # Plot regression fits
  gg.obj <- ggplot() + geom_line(mapping=aes(x=age, y=dnm), data=lines.df[["Pat"]], color=fill.col["Pat"], alpha=1)  
  gg.obj <- gg.obj   + geom_line(mapping=aes(x=age, y=dnm), data=lines.df[["Mat"]], color=fill.col["Mat"], alpha=1)
  gg.obj <- gg.obj + geom_ribbon(mapping=aes(x=age, ymin=lwr, ymax=upr), data=lines.df[["Pat"]], fill=fill.col["Pat"], alpha=0.3)  
  gg.obj <- gg.obj + geom_ribbon(mapping=aes(x=age, ymin=lwr, ymax=upr), data=lines.df[["Mat"]], fill=fill.col["Mat"], alpha=0.3)
  if(!is.null(ref.coef)){ # Include age effects from literature if coef. provided
    for(p in names(lines.df)){
      ref.df <- data.frame("age"=range(lines.df[[p]][["age"]]))
      ref.df[["dnm"]] <- ref.coef[1,p]  +  ref.coef[2,p] * ref.df[["age"]]
      gg.obj <- gg.obj + geom_line(mapping=aes(x=age, y=dnm), data=ref.df, color="#606060", alpha=1, linetype=2)
    }
  }
  
  if(color.f1s){ # Provide each F1 with individual color
    gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm, fill=f1), data=pts.df[["Pat"]], shape=21, size=2.2)
    gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm, fill=f1), data=pts.df[["Mat"]], shape=24, size=2.2) 
  } else {       # Same color for all F1s
    gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm), data=pts.df[["Pat"]], fill=fill.col["Pat"], 
                                  color=stroke.col["Pat"], shape=21, size=2.2)
    gg.obj <- gg.obj + geom_point(aes(x=age, y=dnm), data=pts.df[["Mat"]], fill=fill.col["Mat"], 
                                  color=stroke.col["Mat"], shape=24, size=2.2)
  }
  gg.obj <- gg.obj + xlab("Parental Age At Conception (y)") + 
            ylab("Inferred Number of DNMs") + 
            guides(fill = guide_legend(override.aes=list(shape=21)))
  gg.obj <- gg.obj + ggtitle(use.title) + bw.theme + axis.theme
  return(gg.obj)
}

### Function for printing quantiles
PrintQuantiles <- function(quant, append=TRUE, file=""){
  nm.str <- ""
  q.str <- ""
  for(i in seq_along(quant)){
    nm <- names(quant)[i]
    q  <- quant[i]
    
    nm.str <- sprintf("%s%12s", nm.str, nm)
    q.str  <- sprintf("%s%12.2e", q.str,  q)
  }
  cat(nm.str, "\n", sep="", append=append, file=file)
  cat(q.str, "\n",  sep="", append=append, file=file)
}

##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]

  # Remove non-orthologous DNMs
  dnm.df[[s]] <- dnm.df[[s]][ORTHO == 1]
}


##### LOAD GLM DATA #####
glm.fn <- JoinPath(res.dir, "glm.ortho.RData")
StopNoGLMData(glm.fn)

load(glm.fn) # Loads object save.list
glm.obj <- save.list$glm.obj
glm.sim <- save.list$glm.sim

##### PLOT REGRESSION RESULTS (FIG. 2A-B) #####
ref.coef.ls <- list()
for(s in species){
  # Note that genome size scaling is required to ensure a 1:1 comparison of the slopes
  ref.coef.ls[[s]] <- gao.est * auto.size[s] / jon.size
}

for(s in species){
  cur.title <- paste0(species.nms[s], ": Parental Age Effects")
  
  use.stroke.col <- rep("#333333", 2); names(use.stroke.col) <- parents
  gg.obj <- GGplotAgeEffects(glm.obj[[s]], use.title=cur.title, color.f1s=FALSE, 
                             ref.coef=ref.coef.ls[[s]], fill.col=par.color.vec, stroke.col=use.stroke.col)
  if(s == "hum"){
    exp.plt.df <- data.frame("age"=SeqRange(glm.obj[[s]][["Mat"]]$data$age))
    exp.plt.df[["dnm"]] <- ExponentialModel(exp.plt.df, avg.pub[[s]]["Mat"]) * auto.size[s] / jon.size
    gg.obj <- InsertLayer(gg.obj, after=-3, geom_line(aes(x=age, y=dnm), data=exp.plt.df, color="#606060"))
  }
  
  print(gg.obj)
  
  # Save plot
  out.fn <- JoinPath(res.dir, ifelse(s == "hum", "fig2A.pdf", "fig2B.pdf"))
  if(save.plots){
    cat("Saving", out.fn, "\n")
    ggsave(gg.obj, filename=out.fn, width=5, height=4, units="in", useDingbats=FALSE)
  }
}


##### COMPILE MUTATION TYPES FOR TRANSMITTED DNMS #####
spec.df <- list()
for(s in species){
  # Subset to just transmitted DNMs
  cur.df <- dnm.df[[s]]

  xmit.mask <- (cur.df$TRANSMIT_F2A == 1) | (cur.df$TRANSMIT_F2B == 1)
  xmit.mask[is.na(xmit.mask)] <- F

  cur.df <- cur.df[xmit.mask,]
  
  # Tabulate by mutation type
  spec.vec <- as.character(cur.df$TYPE)
  spec.tab <- table(spec.vec)
    
  spec.df[[s]] <- data.table(spec.tab)
  setnames(spec.df[[s]], old=c("spec.vec", "N"), new=c("TYPE", "COUNT"))
  spec.df[[s]][, PROP := COUNT / sum(COUNT)]
  spec.df[[s]][, Dataset := tolower(species.nms[s])]
  setkey(spec.df[[s]], TYPE)
  
  # Compute confidence intervals
  grp.vec <- cur.df$BLOCK

  prop.ci <- BootCIByGroup(grp.vec, spec.vec)
    
  spec.df[[s]][, LWR := prop.ci[1, TYPE]]
  spec.df[[s]][, UPR := prop.ci[2, TYPE]]
}

# Compute confidence intervals for reference spectra
for(s in species){
  prop.ci <- BootCIMultinom( ref.spec[[s]][["COUNT"]] )
  ref.spec[[s]][, LWR := prop.ci[1,]]
  ref.spec[[s]][, UPR := prop.ci[2,]]
  
  # Add dataset names
  if(s == "hum"){
    ref.spec[[s]][, Dataset := "Jónsson et al."]
  } else if (s == "bab") {
    ref.spec[[s]][, Dataset := "baboon SNPs"]
  }
}

##### TEST MUTATION TYPE PROPORTIONS #####
mut.types <- spec.df$hum[,TYPE]
cur.pval <- FwdVarChiSq(spec.df$hum[mut.types, COUNT], spec.df$bab[mut.types, COUNT], mut.types)
cat("Differences in Human and Baboon spectra (DNM data):\n")
print(cur.pval)

cur.pval <- FwdVarChiSq(ref.spec$hum[mut.types, COUNT], ref.spec$bab[mut.types, COUNT], mut.types)
cat("Differences in Human and Baboon spectra (reference data):\n")
print(cur.pval)


##### PLOT SPECTRA (FIG. 2C) #####
use.colors <- sp.color.vec
names(use.colors) <- tolower(species.nms[species])
use.colors["Jónsson et al."] <- par.color.vec["Pat"]
use.colors["baboon SNPs"]    <- par.color.vec["Mat"]

plt.df <- rbind(spec.df$hum, spec.df$bab,
                ref.spec$hum, ref.spec$bab)

plt.df[, x := match(TYPE, spectra)]
plt.df[, Dataset := factor(Dataset, levels=names(use.colors))]

gg.obj <- ggplot(plt.df) + geom_pointrange(aes(x=x, y=PROP, ymin=LWR, ymax=UPR, color=Dataset), 
                                           shape=21, size=1, fatten=2, fill="white",
                                           position=position_dodge(width=2/3))
gg.obj <- gg.obj + ylab("Proportion of DNMs")
gg.obj <- gg.obj + scale_x_continuous(name="", breaks=seq_along(spectra), labels=spectra)
gg.obj <- gg.obj + scale_color_manual(values=use.colors)

use.title <- c("Mutation Spectrum of Transmitted DNMs\n")

gg.obj <- gg.obj + ggtitle(use.title) + bw.theme + axis.theme  + 
  theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_line(color="grey92",size=0.75),
        legend.background=element_rect(fill="white", size=0.5, linetype="solid", color="grey20"))

gg.obj <- gg.obj + theme(legend.position = c(0.12, 0.75))

print(gg.obj)
# Save plot
out.fn <- JoinPath(res.dir, "fig2C.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=length(spectra)*25/24, height=4, units="in", useDingbats=FALSE)
}

##### LR TEST FOR HIGH PATERNAL AGE EFFECT IN BABOONS #####
if(save.plots){
  lrt.fn <- JoinPath(res.dir, "report.lr_tests.ortho.txt")
  cat("Saving", lrt.fn, "\n")
} else {
  lrt.fn <- ""
}
cat("############## Likelihood Ratio Testing Results ##############\n", append=FALSE, file=lrt.fn)

fit.free <- glm.obj$bab$Pat
dat <- fit.free$data[,c("age","dnm")]

b0.range <- seq(-10, 35, length.out=100) # Range of 
fit.fixd.ls <- list(); fit.pval <- rep(0, length(b0.range))
for(b1.val in c(ref.coef.ls$bab[2,"Pat"], 16/11 * ref.coef.ls$bab[2,"Pat"])){
  for(i in seq_along(b0.range)){
    b0 <- b0.range[i]
    LL.pat.age <- function(b0=b0, b1=1){
      out.ll <- -sum( stats::dpois(dnm, lambda = b0 + b1*age, log=TRUE) )
      return(out.ll)
    }
    fit.fixd.ls[[i]] <- mle2(LL.pat.age, data=dat, fixed=list("b1"=b1.val), method="BFGS") 
    fit.pval[i] <- pchisq(-2*(logLik(fit.fixd.ls[[i]]) - logLik(fit.free)), df=1, lower.tail = FALSE)[1]
  }
  
  cat(sprintf("Using paternal age effect: %0.2f", b1.val), "\n", append=TRUE, file=lrt.fn)
  cat("Quantiles of LR test p-values over range of intercept values\n", append=TRUE, file=lrt.fn)
  PrintQuantiles(quantile(fit.pval, probs=c(0.025, 0.5, 0.975)), append=TRUE, file=lrt.fn)
  
  cat(sprintf("LR test p-value = %0.2e", fit.pval[which(b0.range > 15)[1]]), "\n", append=TRUE, file=lrt.fn)
  
}


##### LR TEST FOR ALPHA #####
cat("\n########################################################\n", append=TRUE, file=lrt.fn)
# alpha at the pat|mat ages of 32.0|28.2 using Gao et al. 3-gen estimates
ref.alpha <- CalcRateAtAge(gao.est[,"Pat"], age=avg.gen$hum["Pat"], denom=1) / 
             CalcRateAtAge(gao.est[,"Mat"], age=avg.gen$hum["Mat"], denom=1)
ref.alpha <- signif(ref.alpha, 4)

cd.ratio <- signif(474 / 217, 2)  # Ratio of cell divisions in human vs baboon sperm
  
a.range <- seq(1,10, length.out=100) # Range of alpha values
fit.pval <- rep(0, length(a.range))
fit.free.ls <- list()

for(i in seq_along(a.range)){
  a <- a.range[i]
  LL.alpha <- function(a = a, pat.b1=1, mat.b0=5, mat.b1=0.5){
    pat.b0 <- a*(mat.b0 + mat.b1*avg.gen$bab["Mat"]) - pat.b1*avg.gen$bab["Pat"]
    pat.lam <- pat.b0 + pat.b1*Pat.age
    mat.lam <- mat.b0 + mat.b1*Mat.age
    
    out.ll <- -sum( stats::dpois(Pat.dnm, lambda = pat.lam, log=TRUE) ) -
      sum( stats::dpois(Mat.dnm, lambda = mat.lam, log=TRUE) )
    return(out.ll)
  }
  
  dat <- list()
  cur.glm <- glm.obj$bab
  for(p in parents){
    for(nm in c("age","dnm")){
      dat[[paste0(p, ".", nm)]] <- cur.glm[[p]]$data[[nm]]
    }
  }
  
  
  fit.free.ls[[i]] <- mle2(LL.alpha, data=dat, method="BFGS") # no constraint on alpha
  if(i == 1){
    fit.fixd <- mle2(LL.alpha, data=dat, fixed=list("a"=ref.alpha / cd.ratio), method="BFGS") 
  }
  
  fit.pval[i] <- pchisq(-2*(logLik(fit.fixd) - logLik(fit.free.ls[[i]])), df=1, lower.tail = FALSE)[1]
}
cat("LR Test of lower alpha in baboons:\n", append=TRUE, file=lrt.fn)
cat("Quantiles of p-values over range of start values\n", append=TRUE, file=lrt.fn)
PrintQuantiles(quantile(fit.pval, probs=c(0.025, 0.5, 0.975)), append=TRUE, file=lrt.fn)


##### TEST FOR SPECIES EFFECT ON SLOPE, INTERCEPT IN EACH SEX ####
# Load puberty regression data
pub.fn <- JoinPath(res.dir, "glm.ortho.postpub.RData")
StopNoGLMData(pub.fn)

load(glm.fn) # Loads object save.list
pub.obj <- save.list$glm.obj
pub.sim <- save.list$glm.sim


LikelihoodRatioTest <- function(h0, h1){
  CalcLikelihood <- function(in.obj){
    if("logLik" %in% class(in.obj)){
      return(in.obj)
    } else if("glm" %in% class(in.obj)){
      return(logLik(in.obj))
    } else {
      stop("Class of in.obj is neither 'logLik' nor 'glm'")
    }
  }
  
  ll.h0 <- CalcLikelihood(h0)
  ll.h1 <- CalcLikelihood(h1)
  
  out.val <- list()
  out.val[["chisq.df"]]   <- attributes(ll.h1)$df - attributes(ll.h0)$df
  out.val[["chisq.stat"]] <- -2*(ll.h0[1] - ll.h1[1])
  out.val[["p.value"]]    <- pchisq(out.val[["chisq.stat"]], df=out.val[["chisq.df"]], lower.tail=FALSE)
  return(out.val)
}
lifehist.timepts <- c("conception", "puberty")
for(lh.tp in lifehist.timepts){
  cat("\n########################################################\n", append=TRUE, file=lrt.fn)
  cat("Testing for species effect on slope (age effect) and intercept in each sex...\n\n", append=TRUE, file=lrt.fn)
  cat("==========================\n", append=TRUE, file=lrt.fn)
  cat(sprintf("Using years since %s\n\n", lh.tp), append=TRUE, file=lrt.fn)
  fix.glm <- list()
  fix.glm[["slope"]] <- list() # Same slope
  fix.glm[["int"]]   <- list() # Same intercept
  fix.glm[["both"]]  <- list() # Same slope and intercept
  for(p in parents){
    sp.df <- data.table()
    for(s in species){
      tmp.df <- as.data.table( glm.obj[[s]][[p]]$data[, c("age", "dnm")] )
      tmp.df[, sp := s]
      sp.df <- rbindlist(list(sp.df, tmp.df))
    }
    if(lh.tp == "puberty"){
      for(s in species){
        sp.df[sp == s, age := age - avg.pub[[s]][p]]
      }
    }
    fix.glm[["slope"]][[p]] <- glm(dnm ~ age + sp,     data=sp.df, family=poisson(link="identity"))
    fix.glm[["int"]][[p]]   <- glm(dnm ~ age:sp + age, data=sp.df, family=poisson(link="identity"))
    fix.glm[["both"]][[p]]  <- glm(dnm ~ age,          data=sp.df, family=poisson(link="identity"))
  }
  
  for(p in parents){
    h1 <- logLik(glm.obj$hum[[p]]) + logLik(glm.obj$bab[[p]])
    attr(h1, "df") <- 4
    for(nm in names(fix.glm)){
      lrt.res <- LikelihoodRatioTest(fix.glm[[nm]][[p]], h1)
      
      out.str <- sprintf("P-value for fixed %5s parameter(s) vs. free model (%3sernal DNMs)", nm, p)
      cat(sprintf("%s: %0.2e\n", out.str, lrt.res$p.value), append=TRUE, file=lrt.fn)
    }
    cat("\n")
  }
}

##### TEST IF MATERNAL MUTATION RATES CONSISTENT WITH RATIO OF AGES ####
cat("\n########################################################\n", append=TRUE, file=lrt.fn)
age.ratio <- unname(avg.gen$hum["Mat"] / avg.gen$bab["Mat"])

age.ratio.range <- seq(0.5, 5, length.out=100) # Range of maternal age ratios to test
fit.pval <- rep(0, length(a.range))
fit.free.ls <- list()

dat <- list()
for(s in species){
  cur.glm <- glm.obj[[s]][["Mat"]]
  for(nm in c("age","dnm")){
    dat[[paste0(s, ".", nm)]] <- cur.glm$data[[nm]]
  }
}

for(i in seq_along(age.ratio.range)){
  r <- age.ratio.range[i]
  LL.mat <- function(r = r, hum.b1=0, bab.b0=1, bab.b1=0){
    hum.b0 <- r*(bab.b0 + bab.b1*avg.gen$bab["Mat"]) - hum.b1*avg.gen$hum["Mat"]
    hum.lam <- hum.b0 + hum.b1*hum.age
    bab.lam <- bab.b0 + bab.b1*bab.age
    
    out.ll <- -sum( stats::dpois(hum.dnm, lambda = hum.lam, log=TRUE) ) -
      sum( stats::dpois(bab.dnm, lambda = bab.lam, log=TRUE) )
    return(out.ll)
  }
  
  
  fit.free.ls[[i]] <- mle2(LL.mat, data=dat, method="BFGS") # no constraint on alpha
  if(i == 1){
    fit.fixd <- mle2(LL.mat, data=dat, fixed=list("r"=age.ratio), method="BFGS") 
  }
  
  fit.pval[i] <- pchisq(-2*(logLik(fit.fixd) - logLik(fit.free.ls[[i]])), df=1, lower.tail = FALSE)[1]
}
cat(sprintf("LR Test of %0.1f-fold faster mutation rate in human mothers:\n", age.ratio), append=TRUE, file=lrt.fn)
cat("Quantiles of p-values over range of start values\n", append=TRUE, file=lrt.fn)
PrintQuantiles(quantile(fit.pval, probs=c(0.025, 0.5, 0.975)), append=TRUE, file=lrt.fn)


##### TEST IF MUTATIONS AT PUBERTY IN FATHERS IS CONSISTENT WITH RATIO OF AGES AT PUBERTY ####
cat("\n########################################################\n", append=TRUE, file=lrt.fn)
pub.ratio <- unname(avg.pub$hum["Pat"] / avg.pub$bab["Pat"])

pub.ratio.range <- seq(1, 5, length.out=100) # Range of paternal age ratios to test
fit.pval <- rep(0, length(a.range))
fit.free.ls <- list()

dat <- list()
for(s in species){
  cur.glm <- pub.obj[[s]][["Pat"]]
  for(nm in c("age","dnm")){
    dat[[paste0(s, ".", nm)]] <- cur.glm$data[[nm]]
  }
}

for(i in seq_along(pub.ratio.range)){
  r <- pub.ratio.range[i]
  LL.pub <- function(r = r, hum.b1=1, bab.b0=1, bab.b1=1){
    hum.b0 <- bab.b0*r
    hum.lam <- hum.b0 + hum.b1*hum.age
    bab.lam <- bab.b0 + bab.b1*bab.age
    
    out.ll <- -sum( stats::dpois(hum.dnm, lambda = hum.lam, log=TRUE) ) -
      sum( stats::dpois(bab.dnm, lambda = bab.lam, log=TRUE) )
    return(out.ll)
  }
  
  
  fit.free.ls[[i]] <- mle2(LL.pub, data=dat, method="BFGS") # no constraint on alpha
  if(i == 1){
    fit.fixd <- mle2(LL.pub, data=dat, fixed=list("r"=pub.ratio), method="BFGS") 
  }
  
  fit.pval[i] <- pchisq(-2*(logLik(fit.fixd) - logLik(fit.free.ls[[i]])), df=1, lower.tail = FALSE)[1]
}
cat(sprintf("LR Test of %0.1f-fold more DNMs at puberty in human fathers:\n", pub.ratio), append=TRUE, file=lrt.fn)
cat("Quantiles of p-values over range of start values\n", append=TRUE, file=lrt.fn)
PrintQuantiles(quantile(fit.pval, probs=c(0.025, 0.5, 0.975)), append=TRUE, file=lrt.fn)
