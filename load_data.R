#!/usr/bin/env Rscript

# Helper script for loading trio and DNM data into R environment

##### READ IN DATA #####
# (as data.table objects)

dat.dir <- "data" # Data directory
res.dir <- "res"  # Results directory

### Data table with analysis modes
mod.df <- LoadModesData("analysis_modes.tsv")

### Data table with trio information
tri.df <- LoadTrioData(JoinPath(dat.dir, "trio.tsv"))

### Data table with DNMs
dnm.df <- list()
map.df <- list()
for(s in species){
  dnm.df[[s]] <- LoadMutData(JoinPath(dat.dir, paste0(s, ".dnms.tsv")))
}

### 50 cM blocks 
map.df <- list()
map.df[["hum"]] <- LoadMapData(JoinPath(dat.dir, "b37.50cM_windows.bed")) # from HapMap b37 
map.df[["bab"]] <- LoadMapData(JoinPath(dat.dir, "Robinson.50cM_windows.bed")) # from Robinson et al

### Reference DNM/SNP Spectra
ref.spec <- list()
ref.spec[["hum"]] <- fread(JoinPath(dat.dir, "hum.dnm_spectra.ortho.tsv"), skip=0)
ref.spec[["bab"]] <- fread(JoinPath(dat.dir, "bab.snp_spectra.ortho.tsv"), skip=0)
for(s in species){
  ref.spec[[s]][, PROP := COUNT/sum(COUNT)]
  setkey(ref.spec[[s]], TYPE)
}

### Table 1 (Alpha across species)
alpha.tab <- fread(JoinPath(dat.dir, "table1.tsv"))
setnames(alpha.tab, 
         old=c("Common Name",
               "Mean Paternal Age at Conception", 
               "Ratio of Paternal to Maternal DNMs"),
         new=c("Common", "Mean_Pat_Age", "Ratio"))
setkey(alpha.tab, Species)

##### PREPROCESS DNM DATAFRAME #####
for(s in species){
  # Add column for genome block group for later bootstrapping
  AddBlockGroup(dnm.df[[s]], s)

  # Identify and mark DNM clusters
  AssignClusters(dnm.df[[s]])
}
