# Day 1 Example Plot

#import the file:
viremia <- read.csv("viremia_data_full.csv")

#view the data in a note:
View(viremia)

#name columns
colnames(viremia) <- c("Bird","n","Family","Order","1","3","4","6")

#choose some colors
cols <-c("black","grey",rainbow(26)[4:26])

#Plot by species - did not work
plot(c(1,3,4,6),as.numberic(viremia[1,6:9]),
  type = "l",lwd =2, ylim = range(viremia [,6:9],na.rm=TRUE),
  xlab = "Day Postinfection",
  ylab = "Log PFU/ml Serum")
for (i in 2:nrow(viremia)){
  lines(c(1,3,4,6),as.numeric(viremia[i,6:9]),lwd =2, col = cols)
}

gc()
------------------------------------------------------------------
# Day 2 Objects and Data Types

# an object stores information in R (an object oriented language)
x <- 5 # numeric
y <- "cat" # character - object y is assigned a character value "cat"
z <- TRUE # logical - only TRUE (T) or FALSE (F) 

# code is always from top to bottom, left to right

# numeric = numbers
# character = text
# logical = TRUE/FALSE 
# factor = fixed categories (ie sex male/female, pets cat/dog/fish/bir)

#data types change how R interprets and uses data (ie "5"+1->error - character can't add to number)

#console: where code runs immediately
#script (.Rmd or .R): where code is written, saved, reused
#best practice: test in console -> save in script/RmD

#Practice:
> a<-5
> b<-"5"
> c<-TRUE
> d<-factor(c("A","B","A","C"))
> a+2
[1] 7
> b == 5 #just checking what was already stated/if a given condition is true
[1] TRUE
> !c #Bang operator (!) means everything but or the opposite
[1] FALSE
> levels(d)
[1] "A" "B" "C"

#"I use AI to help with error correction" - how we are allowed to use AI 
#can't put spaces in object names

# is called commenting

----------------------------------------------------------------------  
#Day 3

#Rmd stands for R mark down
#reproducible reports=text+code+output - a file that saves text
#headers (# in markdown)
#code chunks ('''{r})
#commenting in/out of code

#write code in grey and text outside
# "##" is headers
#outline shows you headers
#Knit "prints" it
  
#output in the top section changes how it is downloaded (PDF, word, HTML)

---------------------------------------------------------------------------

#Day 4
# <- and = are the same, == checks and results in "true"
#for factor, you have to say factor()
# "4" is different than 4 which is different than 4.0 - testing them against each other would result in FALSE

#R Markdown
  #a format for writing reproducible, dynamic reports with R - use it to embed code and results into slideshows, pdfs, etc
  #open, write, embed and render (replace R code with its output and transform the report into whatever you want)
  #menu: File - New File - R Markdown - class of output - OK
  #YAML header is title, author, output, etc.
  #use # for headers: # is big then ## smaller then ### smallest
  #Ctrl+Alt+I puts in 3 ticks (code chunks), subract to 1 for incline code
    #code chunks show the code in the final report
    #inline code replaces the code with the results in the final report
  #commenting outside of code chunk only adds the #, if you want it to be emitted in final report: <!-- your comment -->
  #in the code chunk, the comment will be visible with the #
  #Ctrl+Shift+C is shortcut for commenting
  #eval(or echo)=FALSE - this will leave out the code but will show the output

  
#Day 5
  #syntax
  my_function<-function(parameters) {
    #code to execute.
  }
  my_function(parameters)
  
  odd_or_even <- function(num) {
  }
    
  ?round()
  ?c()
  ?factor()
  ?plot()
  
  #? is help menu, made by random people
  #round
  #x and digits (arguments) are required - goes in ()
  #default is what it decides you want if you don't specify with an argument
  round(pi) #or
  round(x=3.141593) #define x and leave digits default (0)
  round(x=3.141593,digits=3) #digits are after the decimal points
  #homework: go through helps for each of the above and understand them
  
  
  #Day 6 
  #practice for Day 5
  round(4,digits=3)
  round(4.39872, digits=3)
  ?round()
  round(4.39872)
  ?mean  
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,NA,2))
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,2))
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,2),trim=0.1)
  mean(c(3,2,4,5,2,1,6,3,2,2,2,2,NA,2),na.rm=TRUE)

#Day 7
  #read.table(), read.csv() -> import data
  #write.table(), write. csv() -> save/export
  #file paths (relative vs absolute)
  #tip: always check your data after import
    #head() peek at first rows
    #view() open spreadsheet view
    #str() 
  
  
?read.csv #typically use commas as separaters
data <- read.csv(file="viremia_data_full.csv")

head(data) #headers
View(data) #needs to be capitalized
str(data) #structure
dim(data) #column and row dimensions

data <- read.csv(file="bloodmeal_for_BIOL3070.csv")

#Day 8 Mini-report Rubric
#total: 25 points - five grading categories, 5 points each
#1. Structure and Completion
  #All sections present (Abstract, Background, Question/Hypothesis, Methods, Discussion, Conclusion, References)
  #Abstract presents and covers full project (written last)
  #Knits without errors
#2. Background and Question Framing
  #Clear background with at least 1 cited scientific reference
  #Connects to course material (e.g. viremia plot)
  #States study question, hypothesis, and prediction clearly
#3. Methods and Analysis
  #Methods written in prose, not just code
  #At lease 2 analyses/plots with reproducible R code
  #Code commented and organized
#4. Interpretation and Discussion 
  #Interpret each analysis/plot in text
  #Tie results back to hypothesis
  #Note at least 1 limitation or uncertainty
#5. Clarity, References, and Reproducibility
  #Clear, professional writing and formatting 
  #All references cited (scientific + AI if used)
  #Knitted .md file shared via GitHub URL
  #Report looks polished and reproducible

#Day 9 - Various
#commit saves a version of the work
  #records what changed along with a short message
#push uploads the save changes (commits) from RStudio/Posit Cloud to your GitHub repository online
  #how you publish your newest version
# load series and platform data from GEO

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)
install.packages("umap")
library(umap)

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE15852", GSEMatrix =TRUE, AnnotGPL=TRUE)
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
groups <- make.names(c("normal tissue","tumor tissue"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE15852", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE15852", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE15852")




#Day Now - Updated R script
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
force=TRUE
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)
force=TRUE
install.packages("umap")
library(umap)

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE15852", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("10101010101010101010101010101010101010101010101010",
               "101010101010101010101010101010101010")
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
groups <- make.names(c("Breast Tumor Tissue","Normal Breast Tissue"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1], groups[2], sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID"))

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Visualize and quality control test results.
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

