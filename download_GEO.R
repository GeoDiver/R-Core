#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Differential Gene Expression Analysis
# Run           : Rscript download_GEO.R --accession GDS5093 --outrdata ~/Desktop/GDS5093.RData
# ---------------------------------------------------------

#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

suppressMessages(library("argparser"))
suppressMessages(library("GEOquery"))
suppressMessages(library("impute"))
suppressMessages(library("pathview")) # for the id2eg function
# Species database
library("org.Ag.eg.db") # "Anopheles" "aga" "eg"
library("org.At.tair.db") # "Arabidopsis" "ath" "tair"
library("org.Bt.eg.db") # "Bovine" "bta" "eg"
library("org.Ce.eg.db") # "Worm" "cel" "eg"
library("org.Cf.eg.db") # "Canine" "cfa" "eg"
library("org.Dm.eg.db") # "Fly" "dme" "eg"
library("org.Dr.eg.db") # "Zebrafish" "dre" "eg"
library("org.EcK12.eg.db") # "E coli strain K12" "eco" "eg"
library("org.EcSakai.eg.db") # "E coli strain Sakai" "ecs" "eg"
library("org.Gg.eg.db") # "Chicken" "gga" "eg"
library("org.Hs.eg.db") # "Human" "hsa" "eg"
library("org.Mm.eg.db") # "Mouse" "mmu" "eg"
library("org.Mmu.eg.db") # "Rhesus" "mcc" "eg"
library("org.Pf.plasmo.db") # "Malaria" "pfa" "orf"
library("org.Pt.eg.db") # "Chimp" "ptr" "eg"
library("org.Rn.eg.db") # "Rat" "rno" "eg"
library("org.Sc.sgd.db") # "Yeast" "sce" "orf"
library("org.Ss.eg.db") # "Pig" "ssc" "eg"
library("org.Xl.eg.db") # "Xenopus" "xla" "eg"
data(korg)
data(bods)

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

parser <- arg_parser("Input GEO Dataset")
parser <- add_argument(parser, "--accession", default = "GSE51808",
                       help = "Accession Number of the GEO Database")
parser <- add_argument(parser, "--outrdata", default = "GSE51808.RData",
                       help = "Full path to the output rData file")
parser <- add_argument(parser, "--geodbDir", default = ".",
                       help = "Full path to the database directory")
argv   <- parse_args(parser)

#############################################################################
#                          Load Functions                                   #
#############################################################################

# auto-detect if data is log transformed
scalable <- function(X) {
  #  produce sample quantiles corresponding to the given probabilities
  qx <- as.numeric(quantile(X, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  logc <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return (logc)
}

#############################################################################
#                        GEO Input                                          #
#############################################################################

gset <- getGEO(argv$accession, GSEMatrix = TRUE, destdir=argv$geodbDir)

if (grepl('^GDS', argv$accession)) {
  eset           <- GDS2eSet(gset, do.log2 = FALSE)
  gene.names     <- as.character(gset@dataTable@table$IDENTIFIER)
  organism       <- as.character(Meta(gset)$sample_organism)
  gpl            <- getGEO(Meta(gset)$platform, destdir=argv$geodbDir)
  featureData    <- gpl@dataTable@table
} else if (grepl('^GSE', argv$accession)) {
  if (length(gset) > 1) idx <- grep(gset@annotation, attr(gse, "names")) else idx <- 1
  eset           <- gset[[idx]]
  gene.names     <- as.character(eset@featureData@data[, "Gene Symbol"])
  organism       <- as.character(eset@featureData@data[, "Species Scientific Name"][1])
  featureData    <- eset@featureData@data
}

X      <- exprs(eset) # Get Expression Data
pData       <- pData(eset)
rownames(X) <- gene.names

# KNN imputation
if (ncol(X) == 2) {
  X <- X[complete.cases(X), ] # KNN does not work when there are only 2 samples
} else {
  X <- X[rowSums(is.na(X)) != ncol(X), ] # remove rows with missing data
}

# remove all zeros
X <- X[rowSums(X != 0) != 0,]

# Replace missing value with calculated KNN value
tryCatch({
  imputation <- impute.knn(X)
  X          <- imputation$data
}, error=function(e) {
  cat("ERROR: Bad dataset: Unable to run KNN imputation on the dataset.", file=stderr())
  quit(save = "no", status = 1, runLast = FALSE)
})

# If not log transformed, do the log2 transformed
if (scalable(X)) {
  X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
  X <- log2(X)
}

# TODO - what happens if species isn't in the korg/ bods db..
organism.scientific.name <-as.character(korg[which(korg[, "scientific.name"] == organism), "kegg.code"])
organism.common.name <- as.character(bods[which(bods[, "kegg code"] == organism.scientific.name), "species"])

#  Convert Gene Symbols to Entrez IDs (For GAGE Script)
entrez.gene.id <- tryCatch({
  if (c('ENTREZ_GENE_ID') %in% names(featureData)) {
    featureData[, 'ENTREZ_GENE_ID']
  } else {
    package <-as.character(bods[which(bods[, "kegg code"] == organism.scientific.name), "package"])
    # Create two column table containing entrez IDs for geodataset
    entrez.id <- id2eg(ids = gene.names, category = "SYMBOL", pkg.name = package,
                       org = as.character(organism.scientific.name))
    entrez.id[,2]
  }
}, warning = function(warning) {
  cat("# Warning ID2EG: ", file=stderr())
  cat(warning$message, file=stderr())
}, error = function(error) {
  cat("# ERROR ID2EG: ", file=stderr())
  cat(error$message, file=stderr())
  cat("\n", file=stderr())
  cat("## This may be because ID2EG does not support this organism: ", file=stderr())
  cat(paste("(", organism.common.name, ", ", organism.scientific.name, ")\n"), file=stderr())
  return (c('FAILED') )
})

if (! is.na(argv$outrdata)) {
  save(X, pData, gene.names, organism, organism.common.name, organism.scientific.name, entrez.gene.id, file = argv$outrdata)
}
