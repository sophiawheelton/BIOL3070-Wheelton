Final RMarkdown
================
Sophia Wheelton
2025-11-12

- [ABSTRACT](#abstract)
- [BACKGROUND](#background)
- [STUDY QUESTION AND HYPOTHESIS](#study-question-and-hypothesis)
  - [QUESTIONS](#questions)
  - [HYPOTHESIS](#hypothesis)
  - [PREDICTION](#prediction)
- [METHODS](#methods)
  - [FIRST ANALYSIS AND PLOT](#first-analysis-and-plot)
  - [SECOND ANALYSIS AND PLOT](#second-analysis-and-plot)
- [DISCUSSION](#discussion)
  - [INTERPRETATION OF FIRST
    ANALYSIS](#interpretation-of-first-analysis)
  - [INTERPRETATION OF SECOND
    ANALYSIS](#interpretation-of-second-analysis)
- [CONCLUSION](#conclusion)
- [REFERENCES](#references)

# ABSTRACT

# BACKGROUND

Cancer stem cells are tumor cells that undergo self-renewal and
differentiation (Jaggupilli and Elkord, 2012). Past research has linked
these cancer stem cells to the initiation and maintenance of a tumor.
The CD24 surface marker is present on cancer stem cells, and is thus
used to identify and isolate these cancer stem cells.

Research investigating the expression of genes related to the CD24
surface marker (called the CD24 gene in this report) will help us better
understand the causal or correlative role this gene plays in breast
cancer.

# STUDY QUESTION AND HYPOTHESIS

## QUESTIONS

How much does gene expression of CD24 differ between human cancer breast
tissue and normal human breast tissue?

## HYPOTHESIS

The expression of CD24 will differ significantly between tumor breast
tissue and normal breast tissue.

## PREDICTION

We predict that the expression of CD24 will be increased in tumor breast
tissue.

# METHODS

Bong et al. collected total RNAs from individuals with breast cancer and
normal tissues. The total RNAs were processed and hybridized on an
oligonucleotide array chip (specifically the Affymetrix genechip U133A)
for 16 hours. There was a total of 86 sets of gene expression data. Of
these 86 sets, 43 were from individuals with breast tumors and 43 were
from individuals with normal breast tissue. The differences in gene
expression between the tumor and normal groups was compared using paired
t-test (Bong et al., 2010). The Bonferroni multiple test correction was
used to adjust the P-values for the false discovery rate.

The data was processed and analyzed using the GEO2R software. Geo2R uses
a variety of R packages - including DESeq2, GEOquery, and limma (Linear
Models for Microarray Analysis) - to perform differential expression
analysis using series matrix data.

To analyze the data using the GEO2R software, we entered the series
accession number (GSE15852) and defined our sample groups (tumor breast
tissue vs normal breast tissue).

## FIRST ANALYSIS AND PLOT

I am still working on the code for the plot. It is having a hard time
knitting.

``` r
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
```

    ## Bioconductor version '3.18' is out-of-date; the current release version '3.22'
    ##   is available with R version '4.5'; see https://bioconductor.org/install

``` r
BiocManager::install("GEOquery")
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## 'help("repositories", package = "BiocManager")' for details.
    ## Replacement repositories:
    ##     CRAN: http://rspm/default/__linux__/focal/latest

    ## Bioconductor version 3.18 (BiocManager 1.30.26), R 4.3.3 (2024-02-29)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'GEOquery'

    ## Installation paths not writeable, unable to update packages
    ##   path: /opt/R/4.3.3/lib/R/library
    ##   packages:
    ##     boot, class, cluster, codetools, foreign, KernSmooth, lattice, nlme, nnet,
    ##     rpart, spatial, survival

``` r
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## 'help("repositories", package = "BiocManager")' for details.
    ## Replacement repositories:
    ##     CRAN: http://rspm/default/__linux__/focal/latest

    ## Bioconductor version 3.18 (BiocManager 1.30.26), R 4.3.3 (2024-02-29)

    ## Warning: package(s) not installed when version(s) same as or greater than current; use
    ##   `force = TRUE` to re-install: 'limma'

    ## Installation paths not writeable, unable to update packages
    ##   path: /opt/R/4.3.3/lib/R/library
    ##   packages:
    ##     boot, class, cluster, codetools, foreign, KernSmooth, lattice, nlme, nnet,
    ##     rpart, spatial, survival

``` r
library(limma)
```

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
install.packages("umap")
```

    ## Installing package into '/cloud/lib/x86_64-pc-linux-gnu-library/4.3'
    ## (as 'lib' is unspecified)

``` r
library(umap)

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE15852", GSEMatrix =TRUE, AnnotGPL=TRUE)
```

    ## Found 1 file(s)

    ## GSE15852_series_matrix.txt.gz

``` r
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("01010101010101010101010101010101010101010101010101",
               "010101010101010101010101010101010101")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Normal Breast Tissue","Tumor Breast Tissue"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID"))

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
```

    ## [1] "Normal.Breast.Tissue-Tumor.Breast.Tissue"

``` r
ct <- 1  # choose contrast of interest
dT <- topTable(fit2, coef = ct, number = Inf, adjust.method = "BH")

# Define thresholds
logFC_cutoff <- 1
adjP_cutoff <- 0.05

# Assign colors based on thresholds
dT$color <- "grey"  # default
dT$color[dT$logFC >  logFC_cutoff & dT$adj.P.Val < adjP_cutoff] <- "red"    # upregulated
dT$color[dT$logFC < -logFC_cutoff & dT$adj.P.Val < adjP_cutoff] <- "blue"   # downregulated

# Clean up the contrast name (remove periods)
contrast_name <- gsub("\\.", " ", colnames(fit2)[ct])
contrast_name <- gsub("-", " and ", contrast_name)

# Draw volcano plot
with(dT, plot(logFC, -log10(adj.P.Val),
              pch = 20,
              col = color,
              main = paste("GSE15852:", contrast_name),
              xlab = "log2 Fold Change",
              ylab = "-log10 Adjusted P-value"))

# Add legend
legend("bottomleft",
       legend = c("Upregulated", "Downregulated", "Not significant"),
       col = c("red", "blue", "grey"),
       pch = 20,
       bty = "n")
```

![](Final-RMarkdown_files/figure-gfm/pressure-1.png)<!-- -->

## SECOND ANALYSIS AND PLOT

# DISCUSSION

## INTERPRETATION OF FIRST ANALYSIS

## INTERPRETATION OF SECOND ANALYSIS

# CONCLUSION

# REFERENCES

1.  Pau Ni IB, Zakaria Z, Muhammad R, Abdullah N et al. Gene expression
    patterns distinguish breast carcinomas from normal breast tissues:
    the Malaysian context. Pathol Res Pract 2010 Apr 15;206(4):223-8.
    PMID: 20097481
2.  Jaggupilli A, Elkord E. Significance of CD44 and CD24 as cancer stem
    cell markers: an enduring ambiguity. Clin Dev Immunol.
    2012;2012:708036. doi: 10.1155/2012/708036. Epub 2012 May 30. PMID:
    22693526; PMCID: PMC3369436.
3.  OpenAI. (2025). ChatGPT (GPT-5) \[Large language model\]. OpenAI.
    <https://chat.openai.com/>
