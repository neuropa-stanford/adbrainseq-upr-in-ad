## R2.2 — GSVA proper (Bioconductor GSVA), donor = unit of inference.
## Same donor pseudobulk + gene sets as the module-score and ssGSEA analyses; 3-way convergence check.
suppressMessages(library(GSVA))

.base <- "/data/adbrainseq/Stanford U/PERK Human AD brain analysis"
.dirs <- list.dirs(.base, recursive = TRUE)
DATA <- .dirs[grepl("data_extraction$", .dirs)][1]
OUT  <- paste0("/data/adbrainseq/Publication/2024_scRNA seq/",
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

di <- read.csv(file.path(DATA, "donor_info.csv"), stringsAsFactors = FALSE)
di <- di[!duplicated(di$Subject), ]
braak <- setNames(as.integer(di$braaksc), di$Subject)
tang  <- setNames(suppressWarnings(as.numeric(di$tangles)), di$Subject)
grp   <- function(b) ifelse(b <= 2, "low", ifelse(b <= 4, "int", "late"))

setfiles <- c("ER-stress (260)"="ERstress_260_geneset.txt", "PERK (31)"="geneset_PERK.txt",
              "IRE1 (32)"="geneset_IRE1.txt", "ATF6 (74)"="geneset_ATF6.txt",
              "ERAD (75)"="geneset_ERAD.txt",
              "Internal control: mRNA transport (91)"="controlset_MRNA_TRANSPORT.txt")
gsets <- lapply(setfiles, function(f) { x <- readLines(file.path(OUT, f)); x[nchar(trimws(x)) > 0] })
names(gsets) <- names(setfiles)

cells <- c("Ex","In","Ast","Mic","Oli","Opc")
lab <- c(Ex="Excitatory neurons", In="Inhibitory neurons", Ast="Astrocytes",
         Mic="Microglia", Oli="Oligodendrocytes", Opc="OPCs")

cohend <- function(a,b){na<-length(a);nb<-length(b);sp<-sqrt(((na-1)*var(a)+(nb-1)*var(b))/(na+nb-2));(mean(a)-mean(b))/sp}
cifun  <- function(a,b){na<-length(a);nb<-length(b);sp2<-((na-1)*var(a)+(nb-1)*var(b))/(na+nb-2);se<-sqrt(sp2*(1/na+1/nb));tc<-qt(.975,na+nb-2);d<-mean(a)-mean(b);c(d-tc*se,d+tc*se)}

res <- data.frame()
for (ct in cells) {
  tab <- read.csv(file.path(DATA, sprintf("mean_DGE_by_donor_cell_type_%s.csv", ct)),
                  stringsAsFactors = FALSE, check.names = FALSE)
  donors <- tab[["donor"]]
  genecols <- setdiff(colnames(tab)[-(1:2)], c("donor","celltype"))   # mirror the Python column choice
  M <- t(as.matrix(sapply(tab[genecols], function(x) suppressWarnings(as.numeric(x)))))
  rownames(M) <- genecols; colnames(M) <- donors
  M[is.na(M)] <- 0
  param <- gsvaParam(exprData = M, geneSets = gsets, kcdf = "Gaussian", minSize = 5, maxSize = 1000)
  gs <- gsva(param)                                    # gene sets x donors
  g <- sapply(donors, function(d) grp(braak[[d]]))
  for (sn in rownames(gs)) {
    sc <- as.numeric(gs[sn, ]); names(sc) <- donors
    a <- sc[g == "late"]; b <- sc[g == "low"]
    if (length(a) >= 3 && length(b) >= 3) {
      p <- suppressWarnings(wilcox.test(a, b)$p.value); cc <- cifun(a, b)
      res <- rbind(res, data.frame(cell_type=ct, cell_label=lab[[ct]], gene_set=sn,
                 comparison="late vs low", n_late=length(a), n_low=length(b),
                 cohens_d=round(cohend(a,b),3), CI95_low=round(cc[1],3), CI95_high=round(cc[2],3),
                 p_wilcoxon=p, stringsAsFactors=FALSE))
    }
    tv <- tang[donors]; ok <- is.finite(tv) & is.finite(sc)
    if (sum(ok) > 10) {
      cc <- suppressWarnings(cor.test(tv[ok], sc[ok]))
      res <- rbind(res, data.frame(cell_type=ct, cell_label=lab[[ct]], gene_set=sn,
                 comparison="tangles correlation", n_late=sum(ok), n_low=NA,
                 cohens_d=round(as.numeric(cc$estimate),3), CI95_low=NA, CI95_high=NA,
                 p_wilcoxon=cc$p.value, stringsAsFactors=FALSE))
    }
  }
}
res$p_adj_BH <- round(p.adjust(res$p_wilcoxon, "BH"), 4)
write.csv(res, file.path(OUT, "R2Q2_GSVA_donorlevel.csv"), row.names = FALSE)

cat("=== GSVA proper: late (V-VI) vs low (I-II) ===\n")
sub <- res[res$comparison == "late vs low", ]
for (i in seq_len(nrow(sub)))
  cat(sprintf("%-20s %-40s d=%+.2f  p=%.3f  BHq=%.3f\n",
      sub$cell_label[i], substr(sub$gene_set[i],1,38), sub$cohens_d[i], sub$p_wilcoxon[i], sub$p_adj_BH[i]))
cat(sprintf("\nwrote R2Q2_GSVA_donorlevel.csv (%d rows)\n", nrow(res)))
