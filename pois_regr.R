#!/usr/bin/env Rscript
rm(list=ls())

# Script for running Poisson Regressions of age versus DNMs
# under different settings

##### SET UP LIBRARIES AND ENVIRONMENT #####
library("arm")
library("digest")

source("setup_env.R")

##### LOAD DATA #####
source("load_data.R")


##### SETTINGS #####
save.data <- TRUE
run.mode <- 2 # Which comparison to run. See analysis_modes.tsv for run modes

restrict.ortho  <- mod.df[run.mode, restrict.ortho]  # Restrict to the orthologous regions
remove.cpgs     <- mod.df[run.mode, remove.cpgs] # Remove CpGs
use.postpub     <- mod.df[run.mode, use.postpub]  # Use years of parental age post-puberty (rather than post-birth)
use.bamsurg.fnr <- mod.df[run.mode, use.bamsurg.fnr] # Use FNRs estimated using bamsurgeon
midfix          <- mod.df[run.mode, midfix] # output files will be of form res/<prefix>.<midfix>.<suffix>

##### FUNCTIONS #####

### Function for printing a report of calculated rates
ReportRates <- function(glm.list, use.ages, genome.size, n.sim=1000, 
                        ci=0.95, n.signif=2, use.seed=NULL, test.coef=NULL,
                        out.fn="", append=TRUE){
  # Construct format strings with the right number of significant digits
  e.fmt.str <- paste0("%27s: %0.", n.signif, "e (%0.", n.signif, "e, %0.", n.signif, "e)\n")
  f.fmt.str <- paste0("%27s: %0.", n.signif, "f (%0.", n.signif, "f, %0.", n.signif, "f)  %s\n")
  
  CatFunc <- function(in.str){
    cat(in.str, file=out.fn, append=append)
  }
  
  CatFunc(paste("Reporting ", ci*100, "% confidence intervals\n", sep=""))
  
  denom <- 2*genome.size # Denominator to convert to per-bp
  
  # Set seed if given
  if(!is.null(use.seed)){
    set.seed(use.seed)
  }
  
  # 
  sim.obj <- list(); mu.g.pop <- c(); mu.g.pop.sim <- list()
  for(p in names(glm.list)){
    cur.glm <- glm.list[[p]]
    b.vals <- coef(cur.glm) # Coefficients
    b.ci <- confint(cur.glm, level=ci) # CIs on coefficients
    
    CatFunc("=====================================================\n")
    CatFunc(paste(p, "ernal regression", "\n", sep=""))
    
    # Simulate regr coefficients
    sim.obj[[p]] <- sim(cur.glm, n.sims=n.sim)
    
    # Compute point estimate of per generation rate at given population age
    mu.g.pop[p] <- CalcRateAtAge(cur.glm, use.ages[p], denom)
    mu.g.pop.sim[[p]] <- apply(coef(sim.obj[[p]]), 1, function(x) CalcRateAtAge(x, age=use.ages[p], denom=denom))
    cur.ci <- quantile(mu.g.pop.sim[[p]], probs=c(1-ci, 1+ci)/2)
    CatFunc(sprintf(e.fmt.str, sprintf("Mut rate per gen at pop. mean %0.1fy", use.ages[p]), mu.g.pop[p], cur.ci[1], cur.ci[2]))
    
    cur.codes <- SignifCodes(cur.glm)
    for(i in seq_along(b.vals)){
      b <- b.vals[i] # Current coefficient
      nm <- names(b.vals)[i]
      cur.ci <- b.ci[i,]
      CatFunc(sprintf(f.fmt.str, nm, b, cur.ci[1], cur.ci[2], cur.codes[i]))
    }
    
    # Compute LRT p-value against deCODE model
    if(!is.null(test.coef)){
      lambda.vals <- test.coef[1,p] + test.coef[2,p] * cur.glm$data$age
      tmp.ll <- sum(dpois(cur.glm$data$dnm, lambda.vals, log=TRUE))
      
      chi.stat <- -2*as.numeric(tmp.ll - logLik(cur.glm))
      p.val <- pchisq(chi.stat, df=2, lower.tail=FALSE)
      CatFunc("\n")
      CatFunc(sprintf("%27s: %0.2e\n", "LRT p-val vs Gao et al. model", p.val))
      
      # Recompute LRT p-value for human, maternal, removing mothers older than 40
      under.forty <- cur.glm$data$age < 40
      if((p == "Mat") && (!all(under.forty))){
        
        # AIC under Gao et al. exponential model
        exp.lamb <- ExponentialModel(cur.glm$data$age, avg.pub[[s]][p]) * auto.size[s] / jon.size
        exp.loglik <- sum(dpois(cur.glm$data$dnm, lambda=exp.lamb, log=TRUE))
        exp.aic <- 2*(3 - exp.loglik)
        
        # AIC under best fit linear model
        lin.lamb <- (gao.est[1,"Mat"] + gao.est[2,"Mat"] * cur.glm$data$age) * auto.size[s] / jon.size
        lin.loglik <- sum(dpois(cur.glm$data$dnm, lambda=lin.lamb, log=TRUE))
        lin.aic <- 2*(2 - lin.loglik)
        
        recalc.glm <- RegressPoisson(cur.glm$data[under.forty,])
        
        lambda.vals <- test.coef[1,p] + test.coef[2,p] * cur.glm$data$age[under.forty]
        tmp.ll <- sum(dpois(cur.glm$data$dnm[under.forty], lambda.vals, log=TRUE))
        chi.stat <- 2*abs(as.numeric(tmp.ll - logLik(recalc.glm)))
        p.val <- pchisq(chi.stat, df=2, lower.tail = F)
        CatFunc("(recalculated removing mothers >= 40y)\n")
        CatFunc(sprintf("%27s: %0.2e\n", "LRT p-val vs Gao et al. model", p.val))
        CatFunc(sprintf("Exponential AIC: %0.2f\n", exp.aic))
        CatFunc(sprintf("     Linear AIC: %0.2f\n", lin.aic))
        CatFunc(sprintf("----------------------\n"))
        CatFunc(sprintf("      Delta AIC: %0.2f\n", exp.aic - lin.aic))
      }
    }
  }
  
  if(length(glm.list) == 2){
    g.fmt.str <- paste0("Sex-averaged mutation rate per %4s: %0.", n.signif, "e (%0.", n.signif, "e, %0.", n.signif, "e)\n")
    
    mu.g <- sum(mu.g.pop)
    mu.g.ci <- quantile(mu.g.pop.sim[[1]] + mu.g.pop.sim[[2]], probs=c(1-ci, 1+ci)/2)
    
    mu.y <- 2*sum(mu.g.pop) / sum(use.ages)
    mu.y.ci <- quantile(2*(mu.g.pop.sim[[1]] + mu.g.pop.sim[[2]]) / sum(use.ages), probs=c(1-ci, 1+ci)/2)
    
    CatFunc("=====================================================\n")
    CatFunc(sprintf("= Using paternal and maternal ages of %0.1f and %0.1f =\n", use.ages["Pat"], use.ages["Mat"]))
    CatFunc(sprintf(g.fmt.str, "gen",  mu.g, mu.g.ci[1], mu.g.ci[2] ))
    CatFunc(sprintf(g.fmt.str, "year", mu.y, mu.y.ci[1], mu.y.ci[2] ))
  }
  return(invisible(sim.obj))
}


##### PROCESS DATASETS #####
for(s in species){
  # Remove DNMs in clusters of greater than size 2
  dnm.df[[s]] <- dnm.df[[s]][CLUST_SIZE <= 2]

  # Remove non-orthologous DNMs
  if(restrict.ortho){
    dnm.df[[s]] <- dnm.df[[s]][ORTHO == 1]
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

##### INFER DNM COUNTS #####
counts.df <- list()
for(s in species){
  counts.df[[s]] <- list()
  
  d.df <- dnm.df[[s]]
  t.df <- tri.df[s]
  f1.arr <- unique(t.df$F1)
  
  for(p in parents){
    tmp.df <- data.frame()
    for(i in seq_along(f1.arr)){ # Iterate over trios
      f1 <- f1.arr[i]
      f1.row <- t.df[F1 == f1]
      
      cur.age <- f1.row[[paste0(p,"_Age")]][1]
      if( is.na(cur.age) ){  # Skip if no parental age
        cat("Skipping F1", f1, "\n"); next 
      }
      
      # Subtract gestation time from age at birth to obtain ages at conception
      cur.age <- cur.age - avg.con[s]
      if(use.postpub){ 
        # Subtract typical age at puberty to obtain years from puberty >> conception
        cur.age <- cur.age - avg.pub[[s]][p]
      }
      
      tmp.df[f1, "age"] <- cur.age # Set parental age
      
      is.infer <- is.na(f1.row$F2[1])
      
      sub.df <- d.df[F1 == f1]
      if(remove.cpgs){
        d.df <- d.df[TYPE != "CpG>TpG", ]
      }
      
      num.dnm <- nrow(sub.df)
      y1 <- sum(sub.df$CLUST_SIZE == 1)
      y2 <- sum(sub.df$CLUST_SIZE == 2)
      
      if(restrict.ortho){
        use.size <- ortho.size
      } else {
        if(remove.cpgs){
          stop("remove.cpgs requires restrict.ortho to be set to TRUE")
        }
        use.size <- as.numeric(f1.row[, Clb][1])
      }
      
      if(is.infer){ # Transmission info was inferred
        # Proportion read-back phased DNMs that are of the current parent, p
        a <- sub.df[, sum(PHASE_READ == p.chars[p], na.rm=TRUE) / 
                      sum(PHASE_READ %in% p.chars)]
      } else { # Transmission info is known
        if(nrow(f1.row) == 2){ 
          is.xmit <- sub.df[, TRANSMIT_F2A | TRANSMIT_F2B] # Transmitted to either F2
        } else {
          is.xmit <- as.logical(sub.df$TRANSMIT_F2A)
        }
        
        # Number of consensus phased DNMs that are of the current parent, p
        a <- sum((sub.df$PHASE_CONS == p.chars[p]) & is.xmit, na.rm=TRUE) / sum(is.xmit)
      }
      
      if(use.bamsurg.fnr){
        cur.fnr <- mean(tri.df[s, Bamsurg_FNR], na.rm=TRUE)
      } else {
        cur.fnr <- f1.row$FNR[1]
      }
      nu <- a * (y1 * (1-fdr.ones[[s]][f1, est]) + y2 * (1-fdr.twos[[s]])) / 
            ((1 - cur.fnr) * use.size / auto.size[s])
      tmp.df[f1,"dnm"] <- nu
      tmp.df[f1,"F1"]  <- f1
    }
    
    tmp.df <- tmp.df[!is.na(tmp.df$F1),]
    counts.df[[s]][[p]] <- as.data.table(tmp.df)
  }
}

##### RUN POISSON REGRESSION #####
glm.obj <- list() # Store glm objects in a list of lists
for(s in species){
  glm.obj[[s]] <- list()
  for(p in parents){
    glm.obj[[s]][[p]] <- RegressPoisson(counts.df[[s]][[p]])
  }
}

##### PRINT REPORT OF RESULTS #####
ref.coef.ls <- list()
for(s in species){
  # Note that genome size scaling is required to ensure a 1:1 comparison of the slopes
  if(remove.cpgs){
    ref.coef.ls[[s]] <- NULL
  } else {
    ref.coef.ls[[s]] <- gao.est * auto.size[s] / jon.size
  }
}

dig.str <- digest(glm.obj, algo="md5") # Digest of glm.obj
glm.sim <- list()
for(s in species){
# Out file if saving...
  if(save.data){
    rpt.fn <- JoinPath(res.dir, paste("report", midfix, s, "txt", sep="."))
    cat("Saving report to", rpt.fn, "\n")
  } else {
    rpt.fn <- ""
  }

  cat("###########################################################################\n", file=rpt.fn, append=FALSE)
  cat(species.nms[s], "\n", file=rpt.fn, append=TRUE)
  
  # Set seed using first 5 chars of digest string for replicability
  use.seed <- strtoi(paste0("0x", substr(dig.str, 1,5))) + which(species == s) 
  
  test.coef <- ref.coef.ls[[s]]
  use.ages <- avg.gen[[s]]
  if(use.postpub){
    use.ages <- use.ages - avg.pub[[s]] # Shorten reporting ages
    test.coef[1,] <- test.coef[1,] + test.coef[2,]*avg.pub[[s]] # Change intercept
  }
  
  genome.size <- auto.size[s]
  
  glm.sim[[s]] <- ReportRates(glm.obj[[s]], use.ages=use.ages, genome.size=genome.size, use.seed=use.seed,
                              test.coef=test.coef, out.fn=rpt.fn, append=TRUE)  
  
  # Still print to stdout even if saving
  if(rpt.fn != ""){
    writeLines(readLines(rpt.fn))
  }
}

##### SAVE GLM OBJECTS #####

# NOTE: Creating this RData file is necessary for creating some other figures
if(save.data){
  glm.fn <- JoinPath(res.dir, paste0("glm.", midfix, ".RData"))
  
  save.list <- list("glm.obj"=glm.obj, "glm.sim"=glm.sim)
  cat("Saving glm.obj and glm.sim to", glm.fn, "\n")
  save(save.list, file=glm.fn)
}
