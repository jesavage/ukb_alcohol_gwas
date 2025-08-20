# Simulate GWAS results for different interaction model comparisons
# May 26 2024
# Jeanne Savage

# NOTE: Lots of troubleshooting and different programs (GENS2, simuGWAS, simuPOP, sim1000G) were carried out to get this to work. simuGWAS was used in the end based on scripts from Peng & Amos 2010 BMC Bioinfo - https://github.com/BoPeng/simuPOP-examples/blob/master/published/simuGWAS/README.md, but these were substantially modified to get the correct versions of packages and data working and adjust the simulation settings. Use the local versions of all files and a conda environment to ensure replicability.

# set up working environment
cd /simGxE/simuGWAS
conda create --name gens2 python=2.7.18 setuptools numpy=1.10 scipy wxpython
conda activate gens2
conda install -c bpeng simupop

module purge
conda activate gens2

# downloading data, initiating simulation, expanding population, and sampling cases/controls can be run in separate stages by commenting in/out the appropriate lines near the end of the script. this produces a ped file and auxilary phenotype and env/covar files
#python simGxE.py  # 1st sim with large effect sizes, corresponding to output files #1001
python simGxE_2.py 

# create corresponding map file with snp info
awk '{print $2"\t"$1"\t"0"\t"$3}' ex2_init2001.pop.lst > gxe_sample2001.map
sed -i '1d' gxe_sample2001.map

#source scriptSettings.sh
#cd ${TMPDIR}

# prepare files for plink analysis
Rscript - << 'END'
require(data.table)
ped <- fread("gxe_sample2001.ped")
summary(ped[,1:6])
table(ped$V5) #sex
table(ped$V6) #pheno
phen <- fread("gxe_sample2001.pheno")
env <- fread("gxe_sample2001.env")
table(as.numeric(phen$V1))
table(as.numeric(env$V1))
table(as.numeric(phen$V1)+1==ped$V6) #check that phen/env variables stored separately from sim are correctly matched to individuals
cov <- cbind(ped[,1:2],env) #create plink covariate file
names(cov) <- c("FID","IID","env")
fwrite(cov,"gxe_sample2001.cov",sep="\t",quote=F)
eneg <- cov[which(cov$env==0),1:2]
epos <- cov[which(cov$env==1),1:2]
fwrite(eneg,"gxe_sample2001.env_neg",sep="\t",quote=F)
fwrite(epos,"gxe_sample2001.env_pos",sep="\t",quote=F)
END

# check that plink files work and minor alleles of DPL SNPs are correctly identified
plink --file gxe_sample2001 --freq --out test2001

# run plink as in the main analyses for stratified, gxe, and joint 2df tests
## strat env negative
plink --file gxe_sample2001 --keep gxe_sample2001.env_neg --logistic --out gwas2001_strat_envneg
## strat env positive
plink --file gxe_sample2001 --keep gxe_sample2001.env_pos --logistic --out gwas2001_strat_envpos
## gxe individual
plink --file gxe_sample2001 --covar gxe_sample2001.cov --covar-name env --logistic interaction --out gwas2001_gxe
## gxe joing
plink --file gxe_sample2001 --covar gxe_sample2001.cov --covar-name env --logistic interaction --tests 1,3 --out gwas2001_gxe2df

# check results for DPLs
grep rs2044766 *assoc.logistic
grep rs1521011 *assoc.logistic
grep rs7908745 *assoc.logistic


# Manhattan plots
Rscript - << 'END'
require(data.table)
require(qqman)
setwd("/simGxE/simuGWAS/")

envpos <- fread("gwas2001_strat_envpos.assoc.logistic")
envpos <- envpos[is.na(P)==F,]
envneg <- fread("gwas2001_strat_envneg.assoc.logistic")
envneg <- envneg[is.na(P)==F,]
gxe <- fread("gwas2001_gxe.assoc.logistic")
gxem <- gxe[TEST=="ADD" & is.na(P)==F,]
gxei <- gxe[TEST=="ADDxenv" & is.na(P)==F,]
gxe2df <- fread("gwas2001_gxe2df.assoc.logistic")
gxe2df <- gxe2df[TEST=="USER_2DF" & is.na(P)==F,]

envpos$P[envpos$P==0] <- min(envpos$P[envpos$P!=0])
envneg$P[envneg$P==0] <- min(envneg$P[envneg$P!=0])
gxem$P[gxem$P==0] <- min(gxem$P[gxem$P!=0])
gxei$P[gxei$P==0] <- min(gxei$P[gxei$P!=0])
gxe2df$P[gxe2df$P==0] <- min(gxe2df$P[gxe2df$P!=0])


dpl <- c("rs2044766", "rs1521011", "rs7908745")

locmanhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("midnightblue", 
                                                                                   "chartreuse4"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                          annotatePval = NULL, annotateTop = TRUE, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]])) 
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), 
                                    "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
    	
      }
    ticks = c(ticks, d[d$index == i, ]$pos[floor(length(d[d$index ==  i, ]$pos)/2) + 1])
	}
    #ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  # xmax = ceiling(max(d$pos) * 1.03)
  # xmin = floor(max(d$pos) * -0.03)
  xmax = ceiling(max(d$pos) * 1.01)
  xmin = floor(min(d$pos) * .99)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))+15), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "red", pch = 20, 
                             ...))
    text(x=d.highlight$pos, y=d.highlight$logp+8, labels=d.highlight$SNP, cex=.7, col="red")
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos, 
                                                  -log10(P), offset = 0.625, labs = topHits$SNP, 
                                                  cex = 0.45), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos, 
                                                     P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
               labs = topSNPs$SNP, cex = 0.5, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
                  labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}


locmanhattan(envpos, highlight=dpl)
locmanhattan(envneg, highlight=dpl)
locmanhattan(gxem, highlight=dpl)
locmanhattan(gxei, highlight=dpl)
locmanhattan(gxe2df, highlight=dpl)





require(gap)

miamiplot <- function (x, chr = "CHR", bp = "BP", p = "P", pr = "PR", snp = "SNP", 
    col = c("midnightblue", "chartreuse4"), col2 = c("royalblue1", 
        "seagreen1"), ymax = NULL, highlight = NULL, highlight.add = NULL, 
    pch = 19, cex = 0.75, cex.lab = 1, xlab = "Chromosome", ylab = "-log10(P)", 
    lcols = c("black", "red"), lwds = c(3, 1), ltys = c(1, 2), 
    main = "", ...) 
{
    P = index = NULL
    PR = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(pr %in% names(x))) 
        stop(paste("Column", pr, "not found!"))
    if (!(snp %in% names(x))) 
        if (!is.numeric(x[[chr]])) 
            stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    if (!is.numeric(x[[pr]])) 
        stop(paste(pr, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
        PR = x[[pr]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d = subset(d[order(d$CHR, d$BP), ], (P > 0 & P <= 1 & is.numeric(P) & 
        PR > 0 & PR <= 1 & is.numeric(PR)))
    d$logp = -log10(d$P)
    d$logpr = log10(d$PR)
    d$pos = NA
    ymax = ceiling(max(d$logp) + 15)
    ymin = floor(min(d$logpr) - 15)
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }

  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), 
                                    "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
          lastbase
    	
      }
    ticks = c(ticks, d[d$index == i, ]$pos[floor(length(d[d$index ==  i, ]$pos)/2) + 1])
	}
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }

    xmax = ceiling(max(d$pos) * 1.01)
    xmin = floor(min(d$pos) * .99)
    plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        xlim = c(xmin, xmax), ylim = c(ymin, ymax), main = main, 
        xlab = xlab, ylab = ylab, las = 1, pch = pch, cex.lab = cex.lab, 
        ...)
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    col2 = rep(col2, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logpr, cex = cex, pch = pch, ...))
        with(d, points(pos, logp, cex = cex, pch = pch, ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logpr, col = col2[icol], cex = cex, pch = pch, 
                ...))
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], cex = cex, pch = pch, 
                ...))
            icol = icol + 1
        }
    }
    abline(h = -log10(5e-08), col = lcols[2], lwd = lwds[2], 
        lty = ltys[2])
    abline(h = log10(5e-08), col = lcols[2], lwd = lwds[2], lty = ltys[2])
    abline(h = 0, col = lcols[1], lwd = lwds[1], lty = ltys[1])
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "red3", cex = cex, 
            pch = pch, ...))
        with(d.highlight, points(pos, logpr, col = "red3", cex = cex, 
            pch = pch, ...))
    	text(x=d.highlight$pos, y=d.highlight$logp+12, labels=d.highlight$SNP, cex=.7, col="red")
    	text(x=d.highlight$pos, y=d.highlight$logpr-12, labels=d.highlight$SNP, cex=.7, col="red")
    }
    if (!is.null(highlight.add)) {
        print("yessssssssssssssssss")
        if (any(!(highlight.add %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight.add = d[which(d$SNP %in% highlight.add), 
            ]
        with(d.highlight.add, points(pos, logp, col = "darkgreen", 
            cex = cex, pch = pch, ...))
        with(d.highlight.add, points(pos, logpr, col = "darkgreen", 
            cex = cex, pch = pch, ...))
    }
}


envpos$P2 <- envneg$P
miamiplot(x=envpos, pr="P2", highlight=dpl, ylab="-lot10(P) Neg. Subset | -log10(P) Pos. Subset")
gxem$P2 <- gxei$P
miamiplot(x=gxem, pr="P2", highlight=dpl, ylab="-lot10(P) Interaction | -log10(P) Main Effect")

END
