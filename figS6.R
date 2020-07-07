#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plot in Figure SA

##### SET UP LIBRARIES AND ENVIRONMENT #####
# library("arm")
# library("bbmle")
# library("digest")

source("setup_env.R")
set.seed(8080)


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


##### COMPILE MUTATION TYPES IN DNMS #####
spec.df <- list()
for(s in species){
  # Subset to just transmitted DNMs
  cur.df <- dnm.df[[s]]
  
  
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


##### PLOT SPECTRA #####
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

use.title <- c("Mutation Spectrum of All DNMs\n")

gg.obj <- gg.obj + ggtitle(use.title) + bw.theme + axis.theme  + 
  theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_line(color="grey92",size=0.75),
        legend.background=element_rect(fill="white", size=0.5, linetype="solid", color="grey20"))

gg.obj <- gg.obj + theme(legend.position = c(0.12, 0.75))

print(gg.obj)
# Save plot
out.fn <- JoinPath(res.dir, "figS6.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=length(spectra)*25/24, height=4, units="in", useDingbats=FALSE)
}

##### 
mut.types <- spec.df$hum[,TYPE]
cur.pval <- FwdVarChiSq(spec.df$hum[mut.types, COUNT], spec.df$bab[mut.types, COUNT], mut.types)
cat("\nDifferences in Human and Baboon spectra (DNM data):\n")
print(cur.pval)

cur.pval <- FwdVarChiSq(spec.df$hum[mut.types, COUNT], ref.spec$hum[mut.types, COUNT], mut.types)
cat("\nDifferences in Human DNM vs reference spectra:\n")
print(cur.pval)

cur.pval <- FwdVarChiSq(spec.df$bab[mut.types, COUNT], ref.spec$bab[mut.types, COUNT], mut.types)
cat("\nDifferences in Human DNM vs reference spectra:\n")
print(cur.pval)