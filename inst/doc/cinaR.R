## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      cache = FALSE,
                      dev = "png",
                      message = FALSE, error = FALSE, warning = TRUE,
                      fig.dpi = 96)

## -----------------------------------------------------------------------------
library(cinaR)
data("atac_seq_consensus_bm")

## -----------------------------------------------------------------------------
dim(bed)

## -----------------------------------------------------------------------------
# bed formatted file
head(bed[,1:4])

## -----------------------------------------------------------------------------
# create contrast vector which will be compared.
contrasts<- c("B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO", 
              "B6", "B6", "B6", "B6", "B6", "NZO", "NZO", "NZO", "NZO", "NZO", "NZO")

## -----------------------------------------------------------------------------
# If reference genome is not set hg38 will be used!
results <- cinaR(bed, contrasts, reference.genome = "mm10")

## -----------------------------------------------------------------------------
names(results)

## -----------------------------------------------------------------------------
names(results$DA.results)

## -----------------------------------------------------------------------------
colnames(results$DA.results$DA.peaks$B6_NZO)

## -----------------------------------------------------------------------------
head(results$DA.results$DA.peaks$B6_NZO[,1:5])

## -----------------------------------------------------------------------------
head(results$Enrichment.Results$B6_NZO[,c("module.name","overlapping.genes", "adj.p")])

## -----------------------------------------------------------------------------
pca_plot(results, contrasts, show.names = F)

## -----------------------------------------------------------------------------
# Overlaid information
overlaid.info <- c("B6-18mo", "B6-18mo", "B6-18mo", "B6-18mo", "B6-18mo", 
                   "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", "NZO-18mo", 
                   "B6-3mo", "B6-3mo", "B6-3mo", "B6-3mo", "B6-3mo", 
                   "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo", "NZO-3mo")
# Sample IDs
sample.names <- c("S01783", "S01780", "S01781", "S01778", "S01779", 
"S03804", "S03805", "S03806", "S03807", "S03808", 
"S03809", "S04678", "S04679", "S04680", "S04681", 
"S04682", "S10918", "S10916", "S10919", "S10921", 
"S10917", "S10920")

## -----------------------------------------------------------------------------
pca_plot(results, overlaid.info, sample.names)

## -----------------------------------------------------------------------------
heatmap_plot(results)

## -----------------------------------------------------------------------------
heatmap_plot(results, heatmap.peak.count = 200, cluster_cols = F)

## -----------------------------------------------------------------------------
dot_plot(results)

## -----------------------------------------------------------------------------
dot_plot(results, filter.pathways = T)

## -----------------------------------------------------------------------------
contrasts <- sapply(strsplit(colnames(bed), split = "-", fixed = T), 
                    function(x){paste(x[1:4], collapse = "-")})[4:25]
unique(contrasts)

## ----eval=FALSE---------------------------------------------------------------
#  cinaR(..., geneset = new_geneset)

## ----eval=FALSE---------------------------------------------------------------
#  # default geneset to be used
#  data("VP2008")

## ----eval=FALSE---------------------------------------------------------------
#  cinaR(..., batch.correction = T)

## ----eval=FALSE---------------------------------------------------------------
#  # batch information should be number a vector where
#  # the length of it equals to the number of samples.
#  cinaR(..., batch.correction = T, batch.information = c(rep(0, 11), rep(1,11)))

## ----eval=FALSE---------------------------------------------------------------
#  # Ages of the samples could be not in biological interests but should be accounted for!
#  cinaR(..., additional.covariates = c((18, 11), (3, 11)))
#  
#  # More than one covariate for instance, sex and age
#  sex.info <- c("M", "F", "M", "F", "F", "F", "F", "F", "M", "M", "M",
#                "F", "F", "M", "M", "M", "F", "F", "M", "M", "F", "M")
#  
#  age.info <- c((18, 11), (3, 11)
#  covs <- data.frame(Sex = sex.info, Age = age.info)
#  
#  cinaR(..., additional.covariates = covs)

## ----eval=FALSE---------------------------------------------------------------
#  results <- cinaR(..., save.DA.peaks = T, DA.peaks.path = "./Peaks_mice.xlsx")

## ----eval=FALSE---------------------------------------------------------------
#  # new FDR threshold for DA peaks
#  results <- cinaR(..., DA.fdr.threshold = 0.1)
#  
#  # filters out pathways
#  results <- cinaR(..., enrichment.FDR.cutoff = 0.1)
#  
#  # does not run enrichment pipeline
#  results <- cinaR(..., run.enrichment = FALSE)
#  
#  # creates the piechart from chIpSeeker package
#  results <- cinaR(..., show.annotation.pie = TRUE)
#  
#  # change cut-off value for dot plots
#  dot_plot(..., fdr.cutoff = 0.05)

