#!/usr/bin/env Rscript
rm(list=ls())

# Script for generating plots in Figure S4A

##### SET UP LIBRARIES AND ENVIRONMENT #####
source("setup_env.R")
set.seed(888)


##### LOAD DATA #####
source("load_data.R")

res.df <- fread(JoinPath(dat.dir, "tableS1.tsv")) # Sanger results
new.nms <- gsub(names(res.df), pattern=" ", replacement=".")
new.nms[1:4] <- c("CHROM", "POS", "REF", "ALT")
names(res.df) <- new.nms

##### SETTINGS #####
save.plots <- TRUE


##### CONSTANTS #####
s <- "hum"
f1 <- "H14"
ci <- 0.95
quants <- c((1-ci)/2, (1+ci)/2)

##### FIG S4 FUNCTIONS #####


##### PROCESS DATASETS #####
# Restrict dataset to human individual H14
dnm.df <- dnm.df[[s]][F1==f1]
tri.df <- tri.df[s][F1==f1]

# Add F1, cluster, and BLOCK columns to res.df
res.df[, CLUST_ID   := dnm.df$CLUST_ID]
res.df[, CLUST_SIZE := dnm.df$CLUST_SIZE]
res.df[, BLOCK      := dnm.df$BLOCK]

# Remove DNMs in clusters of greater than size 2
dnm.df <- dnm.df[CLUST_SIZE <= 2]
res.df <- res.df[CLUST_SIZE <= 2]


##### SANGER FDR #####
n.dnms <- nrow(res.df)

n.succ <- sum(res.df$Sanger.Result == "Present")
n.fail <- sum(res.df$Sanger.Result == "Absent")
n.sang <- n.succ + n.fail

fdr.sang <- 1 - n.succ / n.sang

sang.ci <- binom.test(x=(n.sang - n.succ), n=n.sang, conf.level=ci)$conf.int


##### TRANSMIT FDR #####
fdr.xmit <- ComputeSingletFDR(tri.df, dnm.df, ci=ci, ignore.clust.size=TRUE)


##### PLOT RESULTS #####
plt.df <- data.frame("x"=1:2, "FDR"=c(fdr.sang, fdr.xmit[,est]),
                     "lwr"=c(sang.ci[1], fdr.xmit[,lwr]),
                     "upr"=c(sang.ci[2], fdr.xmit[,upr]),
                     "Method"=c("Sanger", "Transmission"))
gg.obj <- ggplot(data=plt.df) + geom_pointrange(mapping=aes(x=x, y=FDR, ymin=lwr, ymax=upr, color=Method), shape=21, size=1, fatten=2, fill="white")
gg.obj <- gg.obj + ylim(0,0.5) + scale_x_continuous(name="", breaks=plt.df$x, labels=plt.df$Method, limits=c(0.5,2.5)) + ylab("False Discovery Rate")
gg.obj <- gg.obj + scale_color_manual(values=c("mediumorchid","green4"))
gg.obj <- gg.obj + ggtitle("Comparison of Sanger Seq results with\nFDR estimation from transmission")
gg.obj <- gg.obj + bw.theme + axis.theme + xgrid.off.theme
print(gg.obj)

out.fn <- JoinPath(res.dir, "figS4A.pdf")
if(save.plots){
  cat("Saving", out.fn, "\n")
  ggsave(gg.obj, filename=out.fn, width=4, height=4, units="in", useDingbats=FALSE)
}


##### STATISTICAL TESTS #####
n.xmit <- sum(dnm.df$TRANSMIT_F2A)
two.way.tab <- rbind( c(n.succ, n.fail),
                      c(n.xmit, n.xmit - ceiling(tri.df[,Xmit_Prob]*n.dnms)) )

cat("-------------------------", "\n")
cat("| Sang Succ | Sang Fail |", "\n")
cat("-------------------------", "\n")
cat("| Xmit Obs  | Xmit Miss |", "\n")
cat("-------------------------", "\n")

cat(sprintf("%3d", two.way.tab[1,]), "\n")
cat(sprintf("%3d", two.way.tab[2,]), "\n")

fsh.t <- fisher.test(two.way.tab)

cat("\n-------------------------\n")
cat(sprintf("Fisher's Exact Test p-val: %0.4f\n", fsh.t$p.value))
