#!/usr/bin/env Rscript

# Helper script for setting up R environment with libraries, global constants

##### LIBRARIES #####
library("MASS")
library("data.table")
library("bit64")
library("ggplot2")


##### CONSTANTS #####
species.nms  <- c("Human", "Baboon") # Species names
species <- tolower( substr(x=species.nms, 1, 3) ) # Species short names
names(species.nms) <- species

parents <- c("Pat","Mat") # Strings denoting parents
p.chars <- c("P","M")     # Chars denoting parents
names(p.chars) <- parents

spectra <- c("C>A", "C>G",
             "T>A","T>C","T>G",
             "CpG>TpG","CpH>TpH") # spectra names

##### LIFE HISTORY CONSTANTS #####

# Typical ages at conception (in years)
avg.con <- c(280, 171) / 365
names(avg.con) <- species

# Typical generation time
avg.gen <- list("hum"=c(32.0, 28.2),
                "bab"=c(10.7,10.2))

# Typical age at puberty
avg.pub <- list("hum"=c(13,12),
                "bab"=c(5.41,4.5))

for(s in species){
  names(avg.gen[[s]]) <- parents
  names(avg.pub[[s]]) <- parents
}

##### GENOME CONSTANTS #####

# Size of the (haploid) orthologous genome
ortho.size <- 1631476416

# Size of the (haploid) autosomal genome
auto.size <- c(2881033286, 2581196250)
names(auto.size) <- species

# Size of the (haploid) non-cpg autosomal genome
noncpg.size <- c(2854280589, 2554556385)
names(noncpg.size) <- species

# Size of the (haploid) genome reported by Jonsson et al., Nature 2017
jon.size <- 2682890000

# Size of RepeatMask regions of the genomes
rep.size <- list()
rep.size[["ortho"]] <- c(730854298,  483152809)
rep.size[["auto"]]  <- c(1339775416, 915593622)
for(nm in names(rep.size)){
  names(rep.size[[nm]]) <- species
}

##### REFERENCE MODEL PARAMETERS #####

# Three-generation estimates from Gao, et al. PNAS 2019
gao.est <- cbind(c(7.95, 1.40),
                 c(3.54, 0.33))
colnames(gao.est) <- parents

# Gao et al. estimates for exponential maternal age effect
gao.mat.exp <- c()
gao.mat.exp["a.m"] <- 7.8
gao.mat.exp["b.m"] <- 0.072
gao.mat.exp["c.m"] <- 0.45


##### GGPLOT THEMES, COLORS, FORMATTING #####
sp.color.vec <- c("#00BFC4", "#FFAB00") # Species colors (teal, orange)
names(sp.color.vec) <- species

par.color.vec <- c("#1C7FDB", "#FF4E15") # Parent colors (blue, red)
names(par.color.vec) <- parents


axis.theme <- theme(axis.text=element_text(size=11))
bw.theme <-  theme_bw(base_size = 12, base_family = "ArialMT") + 
             theme(panel.border=element_rect(size = 1.5))
xgrid.off.theme  <- theme(panel.grid.major.x=element_blank(), 
                          panel.grid.minor.x=element_blank())

##### FUNCTIONS #####

### Function for joining paths together
JoinPath <- function(...){
  path.parts <- list(...)

  out.path <- path.parts[[1]] # Out path
  path.parts[[1]] <- NULL
  for(x in path.parts){
    out.path <- paste(out.path, x, sep="/")
  }

  while( grepl(out.path, pattern="//", fixed=TRUE) ){
    out.path <- gsub(x=out.path, pattern="//", replacement="/", fixed=TRUE)
  }

  return(out.path)
}

### Function for loading analysis modes
LoadModesData <- function(fn){
  out.df <- fread(fn, header=TRUE)
  for(nm in names(out.df)){
    if(nm != "midfix"){
      out.df[, eval(nm) := as.logical(get(nm))]
    }
  }
  return(out.df)
}


### Function for loading and preprocessing DNM dataset
LoadMutData <- function(fn){
  col.cls <- c("CHROM"="character", "F1"="character", "F2A"="character", "F2B"="character",
               "PHASE READ"="character", "PHASE PEDIGREE F2A"="character", "PHASE PEDIGREE F2B"="character",
               "PHASE CONSENSUS"="character")
  out.df <- fread(fn, colClasses=col.cls, na.strings=c("."))
  
  # Rename some columns to make code less verbose
  setnames(out.df, old=c("TRANSMIT F2A", "TRANSMIT F2B", 
                         "PHASE READ", "PHASE PEDIGREE F2A", 
                         "PHASE PEDIGREE F2B", "PHASE CONSENSUS",
                         "IN ORTHOLOGOUS", "CPG STATUS", "IN REPEAT"),
                   new=c("TRANSMIT_F2A", "TRANSMIT_F2B", 
                         "PHASE_READ", "PHASE_PED_F2A",
                         "PHASE_PED_F2B", "PHASE_CONS",
                         "ORTHO", "CPG", "REPEAT"))
  setkey(out.df, F1)
  return(out.df)
}

### Function for loading and preprocessing dataset w/ trio-specific info
LoadTrioData <- function(fn){
  col.cls <- c("Callable Genome Size"="numeric", "Callable Repeat Genome Size"="numeric")
  out.df <- fread(fn, na.strings=c("."), colClasses=col.cls)
  setnames(out.df, old=c("Father", "Mother",
                         "Father Age At Birth", "Mother Age At Birth", 
                         "F1 Depth of Coverage", "Bamsurgeon FNR", "Transmission Probability (q tilde)", 
                         "Callable Genome Size", "Callable Repeat Genome Size"),
                   new=c("Pat", "Mat",
                         "Pat_Age", "Mat_Age", 
                         "F1_Depth", "Bamsurg_FNR", "Xmit_Prob", 
                         "Clb", "Clb_Rep"))
  
  ### Process trio dataframe
  out.df[, Notes := NULL] # Remove Notes column
  out.df[, Sps := substr(Species, 1, 3)] # Add shorthand species column 
  setkey(out.df, Sps) # Set "Sps" as key
  
  return(out.df)
}

### Function for loading and preprocessing dataset w/ trio-specific info
LoadMapData <- function(fn){
  out.df <- fread(fn, col.names=c("CHROM","BEGIN","END"),
                  colClasses=c("character", "integer", "integer"))
  out.df[, NAMES := paste(CHROM, BEGIN, sep="_")]
  return(out.df)
}

### Function for adding a column for genome block group to a DNM dataframe
### (for later bootstrapping ), requires species (s)
AddBlockGroup <- function(df, s){
  if("BLOCK" %in% names(df)){
    warning("BLOCK column already present. Skipping AddBlockGroup...")
    return(NULL)
  }
  
  if(!exists("map.df")){
    stop("Cannot AddBlockGroup to DNM data.table without map.df object")
  }
  new.col <- rep("", nrow(df))
  for(i in 1:nrow(df)){
    cur.row <- df[i,]
    new.col[i] <- map.df[[s]][(CHROM == cur.row$CHROM) & 
                              (BEGIN <  cur.row$POS) & 
                              (END   >= cur.row$POS), NAMES]
  }
  df[, BLOCK := new.col]
  
  return(NULL)
}


### Function for assigning clusters to a DNM dataframe
AssignClusters <- function(df, max.dist=100){
  id.col <-    "CLUST_ID"
  count.col <- "CLUST_SIZE"
  if( all(c(id.col, count.col) %in% names(df)) ){ 
    return(NULL) 
  }
  
  MakeClustId <- function(f1, chr, pos){
    return( paste("clust", f1, chr, pos, sep="_") )
  }
  
  ClustPositions <- function(pos, max.dist){ 
    # Clusters a vector of (sorted positions) from a single chromosome
    # Returns a vector denoting the cluster assignment using the first position # in the cluster
    if(length(pos) == 1){ return( pos ) }
    
    out.pos <- rep(0, length(pos))
    out.pos[1] <- pos[1]
    for(i in 2:length(pos)){
      if(pos[i] - pos[i-1] <= max.dist){
        out.pos[i] <- out.pos[i-1]
      } else {
        out.pos[i] <- pos[i]
      }
    }
    return(as.integer(out.pos))
  }
  
  # Assign cluster names
  df[, eval(id.col) := MakeClustId(F1, CHROM, ClustPositions(POS, max.dist)), by=.(F1, CHROM)]
  
  df[, eval(count.col) := rep(length(POS), length(POS)), by=.(CLUST_ID)]
  
  return(NULL)
}

### Function for creating a sequence from min to max
SeqRange <- function(x, length.out=100){
  return( seq(from=min(x), to=max(x), length.out=length.out) )
}

### Insert a new layer into a ggplot
InsertLayer <- function(plt.obj, after=0, ...){
  # plt.obj : Plot object
  # after   : Position where to insert new layers, relative to existing layers
  # ...     : additional layers, separated by commas (,) instead of plus sign (+)

  if(after < 0){
    after <- after + length(plt.obj$layers)
  }

  if(!length(plt.obj$layers)){
    plt.obj$layers <- list(...)
  }
  else{
    plt.obj$layers <- append(plt.obj$layers, list(...), after)
  }

  return(plt.obj)
}


### Function for bootstrap resampling by some group
# (e.g., 50 cM blocks) 
ResampleByGroup <- function(grp.arr, itm.arr, n.sim=1000, grp.names=NULL){
  n <- length(grp.arr)
  if(length(itm.arr) != n){
    stop("grp.arr and itm.arr must have same length")
  }

  itm.uniq <- unique(itm.arr)
  n.itm    <- length(itm.uniq)

  grp.uniq <- unique(grp.arr)
  n.grp    <- length(grp.uniq)

  # Build matrix with rows as groups, columns as items
  mat <- matrix(0, nrow=n.grp, ncol=n.itm)
  rownames(mat) <- grp.uniq
  colnames(mat) <- itm.uniq
  for(g in grp.uniq){
    itm.vals <- table(itm.arr[(grp.arr == g)])
    for(i in names(itm.vals)){
      mat[g,i] <- itm.vals[i]
    }
  }

  if(is.null(grp.names)){
    grp.names <- grp.uniq
  } else {
    grp.names <- unique(grp.names)
    if(!all(grp.uniq %in% grp.names)){
      stop("all values in grp.arr must be in grp.names")
    }

    n.extra <- length(grp.names) - n.grp
    if(n.extra > 0){
      mat.extra <- matrix(0, nrow=n.extra, ncol=n.itm)
      rownames(mat.extra) <- grp.names[!(grp.names %in% grp.uniq)]
      colnames(mat.extra) <- itm.uniq

      mat <- rbind(mat, mat.extra)
    }
  }
  grp.names <- rownames(mat)
  n.row <- nrow(mat)


  sim.vals <- replicate( colSums(mat[sample(grp.names, size=n.row, replace=TRUE),]), n=n.sim )
  return(sim.vals)
}

### Compute confidence intervals by resampling by some group
BootCIByGroup <- function(grp.arr, itm.arr, n.sim=1000, ci=c(0.025,0.975), grp.names=NULL){
  sim.vals <- ResampleByGroup(grp.arr, itm.arr, n.sim, grp.names)
  
  # sim.vals <- replicate( colSums(mat[sample(grp.uniq, size=n.grp, replace=TRUE),]), n=n.sim )
  sim.vals <- apply(sim.vals, 2, function(x) x/sum(x))
  
  q.out <- apply(sim.vals, 1, function(x) quantile(x, probs=ci))
  
  return(q.out)
}

### Compute confidence intervals by multinomial resampling ###
BootCIMultinom <- function(count.vec, n.sim=1000, ci=c(0.025,0.975)){
  n <- sum(count.vec)
  prop <- count.vec / sum(count.vec)

  sim.prop <- rmultinom(n.sim, size=n, prob=prop) / n
  out.tab <- apply(sim.prop, 1, function(x) quantile(x, ci))

  return(out.tab)
}

### Forward variable selection (a.k.a. ordered Chi Square tests)
FwdVarChiSq <- function(pop.a, pop.b, mut.types){
  # arguments should be flat arrays
  if(length(pop.a) != length(pop.b)) { stop("Lengths of all arguments must match") }
  if(length(pop.a) != length(mut.types)) { stop("Lengths of all arguments must match") }
  
  cur.df <- data.table(TYPE  = mut.types,
                       A     = pop.a,
                       B     = pop.b,
                       P_VAL = rep(0.0, length(mut.types)),
                       P_ORD = rep(0.0, length(mut.types)))
  
  # Calculate unordered p-values
  for(i in 1:cur.df[,.N]){
    mat <- cbind(c(cur.df[i,A], sum( cur.df[-i,A] )),
                 c(cur.df[i,B], sum( cur.df[-i,B] )))
    
    cur.df[i, P_VAL := chisq.test(mat)$p.value]
  }
  
  # Calculate ordered p-values
  cur.df <- cur.df[order(P_VAL),]
  
  
  cur.df[1, P_ORD := P_VAL] # Type 1
  for(i in 2:(cur.df[,.N] - 1)){ # Types 2 through (N-1)
    mat <- cbind(c(cur.df[i,A], sum( cur.df[(i+1):.N,A] )),
                 c(cur.df[i,B], sum( cur.df[(i+1):.N,B] )))
    cur.df[i, P_ORD := chisq.test(mat)$p.value]
  }
  
  i <- cur.df[,.N] # Type N
  mat <- cbind(c(cur.df[i,A], sum( cur.df[i-1,A] )),
               c(cur.df[i,B], sum( cur.df[i-1,B] )))
  cur.df[i, P_ORD := chisq.test(mat)$p.value]
  
  setkey(cur.df, TYPE)
  out.p <- cur.df[mut.types, P_ORD]
  names(out.p) <- mut.types
  return(out.p)
}

### Function for obtaining p-values of coefficients from a linear model
ExtractPval <- function(in.model){
  return(summary(in.model)$coefficients[,4])
}

### Function for producing "significance" strings for coefficients from a (general) linear model
SignifCodes <- function(in.model){
  p.vals <- ExtractPval(in.model)
  out.vals <- c()
  for(p in p.vals){
    if(p < 0.001){
      cur.char <- '***'
    } else if(p < 0.01){
      cur.char <- '**'
    } else if(p < 0.05){
      cur.char <- '*'
    } else if(p < 0.1){
      cur.char <- '.'
    } else {
      cur.char <- ' '
    }
       
    out.vals <- c(out.vals, cur.char)
  }
  return(out.vals)
}

### Determine the proper value of 'q', the expected 
# transmission probability for a given F2
ComputeQ <- function(pat.nm, xmit.prob){
  if(pat.nm == "1X2816"){
    q <- 0.5 * xmit.prob / 0.625
  } else {
    q <- xmit.prob
  }
  return(q)
}

### Calculate FDR given transmission observations (z),
# and expected transmission probability (q)
CalcFDR <- function(z, q){
  z.len <- length(z)
  y <- sum(z) # Number of DNMs observed

  if(z.len == 2) {
    if( any(names(z) != c("1", "0")) ){
      stop("z must be ordered: ('1', '0')")
    }

    # Minimum # of true DNMs possible (aka the # of DNMs that were transmitted at least once)
    min.w <- y - z["0"]
    w.vals <- min.w:y # Possible values of w, true number of DNMs
    lhood <- dbinom(x=z["1"], size=w.vals, prob=q, log=TRUE)
  } else if(z.len == 4) {
    if( any(names(z) != c("11", "10", "01", "00")) ){
      stop("z must be ordered: ('11', '10', '01', '00')")
    }

    min.w <- y - z["00"]
    w.vals <- min.w:y # Possible values of w, true number of DNMs

    q.vec <- as.vector( t( c(q[1], 1-q[1]) %o% c(q[2], 1-q[2]) ) ) # Vector of length 4
    w.vals <- min.w:y
    lhood <- rep(NA, length(w.vals))
    for(i in seq_along(w.vals)){
      cur.w <- w.vals[i]
      lhood[i] <- dmultinom(c(z[1:3], cur.w-min.w), size=cur.w, prob=q.vec, log=TRUE)
    }
  } else {
    stop("z must be of length 2 (binomial) or 4 (multinomial)")
  }
  
  out.fdr <- 1 - w.vals[which.max(lhood)] / y

  return(out.fdr)
}

### Compute FDR from a DNM dataframe for SINGLETS
ComputeSingletFDR <- function(t.df, d.df, ci=NULL, ignore.clust.size=FALSE){
  if(!ignore.clust.size){
    d.df <- d.df[CLUST_SIZE == 1]
  }
  
  f1.arr <- unique(t.df$F1)
  fdr.ones <- data.table( "F1"=f1.arr, "est"=as.double(rep(NA, length(f1.arr))) )
  if(!is.null(ci)){
    fdr.ones[, lwr := est]
    fdr.ones[, upr := est]
  }
  setkey(fdr.ones, F1)

  for(f1 in f1.arr){ # For each F1
    cur.row <- t.df[F1 == f1]
    f2.vec <- cur.row[, F2]
    
    # Skip if no F2s
    if( any(is.na(f2.vec)) ){ next }
    
    q <- cur.row[, ComputeQ(Pat, Xmit_Prob), by=F2]$V1

    if(length(f2.vec) == 2){ # If two F2s
      z.strs <- paste(d.df[F1==f1, TRANSMIT_F2A], d.df[F1==f1, TRANSMIT_F2B], sep="")
      z.nms <- c("11", "10", "01", "00")
    } else if(length(f2.vec) == 1){ # If one F2
      z.strs <- as.character(d.df[F1==f1, TRANSMIT_F2A])
      z.nms <- c("1", "0")
    }
    
    z <- table(z.strs)[z.nms]
    fdr.ones[f1, est := CalcFDR(z,q)]

    if(!is.null(ci)){ # If confidence intervals requested, return as well
      z.sim <- ResampleByGroup(d.df[F1==f1, BLOCK], z.strs)
      z.sim <- z.sim[z.nms,]

      tmp.fdr.sim <- apply(z.sim, 2, function(x){ CalcFDR(x,q) })
      fdr.ones[f1, c("lwr", "upr") := 
                   as.list(quantile(tmp.fdr.sim, c(1-ci,1+ci)/2))]
    }
  }
  return(fdr.ones)
}

### Function for computing FDR from a DNM dataframe
# for DOUBLETS
ComputeDoubletFDR <- function(t.df, d.df, ci=NULL){
  d.df <- d.df[CLUST_SIZE == 2 & !is.na(F2A)]
  z <- c(0, 0); names(z) <- c("1", "0")

  for(clust.id in unique(d.df$CLUST_ID)){
    tmp.df <- d.df[CLUST_ID == clust.id]

    xmit.status <- tmp.df[, TRANSMIT_F2A == 1]
    if( !is.na(tmp.df[1,F2B]) ){
      xmit.status <- xmit.status | (tmp.df[, TRANSMIT_F2A == 1])
    }
    z <- z + c(all(xmit.status), !all(xmit.status))
  }

  mean.q <- t.df[!duplicated(F1) & !is.na(F2), ComputeQ(Pat, Xmit_Prob), by=F1]$V1
  mean.q <- mean(mean.q)

  fdr.twos <- CalcFDR(z, mean.q)

  return(fdr.twos)
}

### Function for computing FDR for F1s without F2
RegressFDR <- function(t.df, fdr.df){
  dep.df <- t.df[, .(DEP = F1_Depth[1]), by=F1] # Depth values
  setkey(dep.df, F1)
  fdr.df[, pred  := is.na(est)]      # Add column indicating which FDR est predicted
  fdr.df[, depth := dep.df[F1, DEP]] # Add column indicating FDR depth
  
  # Run linear regression
  fdr.lm <- lm(est ~ depth, data=fdr.df[!(pred)])

  # Predict FDR for F1s without F2 using regression
  fdr.df[(pred), est := predict(fdr.lm, newdata=fdr.df[(pred)])]

  return( invisible(fdr.lm) )
}

### Calculates mutation rate from coefficients or from a GLM object 
### at a given age of reproduction
CalcRateAtAge <- function(in.obj, age, denom){
  if("lm" %in% class(in.obj)){
    # GLM/LM
    out.val <- predict(in.obj, newdata=data.frame("age"=age))
  } else if("numeric" %in% class(in.obj)) { 
    # Coefficients
    if(length(in.obj) != 2){
      stop("Length of 'in.obj' must be 2 if passing single set of coefficients")
    }
    out.val <- unname(in.obj[1] + in.obj[2]*age)
  } else if ("matrix" %in% class(in.obj)){
    # Matrix of coefficients
    if(ncol(in.obj) != 2){
      stop("Number of columns of 'in.obj' must be 2 if passing multiple sets of coefficients")
    }
    out.val <- in.obj[,1] + in.obj[,2]*age
  }
  return(out.val / denom)
}

### Perform Poisson regression
RegressPoisson <- function(df, do.round=TRUE){
  # Names of df should include "dnm", "age"
  df <- data.table::copy(df)
  if(do.round){
    df[, dnm := round(dnm)]
  }
  
  m <- glm(dnm ~ age, data=df, family=poisson(link="identity"))
  return(m)
}

### Calculate the Poisson mean using the Gao et al.
### exponential model parameter estimates
ExponentialModel <- function(age, pub.age, param=gao.mat.exp){
  # age: age at conception
  # pub.age: typical age at puberty
  # Returns lambda parameter (Poisson mean)
  lamb <- param["a.m"] + 
    exp( (age - pub.age)*param["b.m"] + param["c.m"] )
  return(lamb)
}

### Build a tree-formatting string to pass to read.tree()
BuildTreeStr <- function(col, tb) { 
  tree.str <- paste0("((((((human:",tb[1,col],",chimpanzee:",tb[2,col],"):",tb[3,col],",orangutan:",
                     tb[4,col],"):",tb[5,col],",(((rhesus_macaque:",tb[6,col],",crab-eating_macaque:",
                     tb[7,col],"):",tb[8,col],",baboon:",tb[9,col],"):",tb[10,col],",green_monkey:",
                     tb[11,col],"):",tb[12,col],"):",tb[13,col],",(marmoset:",tb[14,col],",squirrel_monkey:",
                     tb[15,col],"):",tb[16,col],"):",tb[17,col],",bushbaby:",tb[18,col],"):",tb[19,col],
                     ",mouse:",tb[20,col],");")
  return(tree.str)
}

### Extract Moorjani et al. substitution rates for given column name
ComputeSubRates <- function(data.path, col.nm, kp.species=c("human", "baboon")){
  # Load substitution data for(average rates)
  tb <- read.table(data.path, header=T, sep="\t")
  
  # Build tree
  col <- which(colnames(tb) == col.nm)
  full.tree <- read.tree( text=BuildTreeStr(col, tb) )
  
  # Remove species that aren't baboon or human
  prun.tree <- keep.tip(full.tree, kp.species)
  
  # Get substitution rates
  out.rates <- distRoot(prun.tree)
  names(out.rates) <- substr(names(out.rates), 1,3)
  return( list('rates'=out.rates, 'tree'=prun.tree) )
}

### Function for computing distance to MRCA for two tip labels of a tree
Dist2MRCA <- function(tree.obj, tip.nms){
  d.root <- c() # Tip to root distances
  for(nm in tip.nms){
    d.root <- c(d.root, distRoot(tree.obj, nm))
  }
  d.tips <- as.matrix(distTips(tree.obj, tips=tip.nms))[tip.nms[1], tip.nms[2]] # Distance between tips
  d.node <- (sum(d.root) - d.tips)/2 # Distance from MRCA node to root
  
  out.dist <- d.root - d.node
  names(out.dist) <- tip.nms
  return(out.dist)
}


### Stop and print a message if GLM Data is missing
StopNoGLMData <- function(glm.fn){
  if(!file.exists(glm.fn)){
    # Which run mode would be appropriate?
    midfix <- gsub(glm.fn, pattern=".*glm\\.", replacement="")
    midfix <- gsub(midfix, pattern="\\.RData", replacement="")
    run.mode <- which(mod.df$midfix %in% midfix)
    
    if(length(run.mode) > 0){
      stop.msg <- sprintf("Can't find file '%s'\nTry running pois_regr.R with the settings:\n", glm.fn)
      stop.msg <- paste0(stop.msg, sprintf("  run.mode <- %d\n", run.mode))
    } else {
      stop.msg <- sprintf("Can't find file '%s'\nUnrecognized file name. Check 'analysis_modes.tsv' for run modes.\n", glm.fn)
    }
    stop(stop.msg)
  }
  cat("GLM data found.\n")
}

##### SET SEED #####
set.seed(8080)
