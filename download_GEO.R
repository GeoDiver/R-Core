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

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

parser <- arg_parser("Input GEO Dataset")
parser <- add_argument(parser, "--geodbpath",
                       help = "GEO Dataset full path")
parser <- add_argument(parser, "--accession", default = "GDS5093",
                       help = "Accession Number of the GEO Database")
parser <- add_argument(parser, "--outrdata",
                       help = "Full path to the output rData file")
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

if (is.na(argv$geodbpath)) {
  gse <- getGEO(argv$accession, GSEMatrix = TRUE)
} else {
  gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE)
}
eset        <- GDS2eSet(gse, do.log2 = FALSE)
pData       <- pData(eset)

X           <- exprs(eset)  # Get Expression Data
gene.names  <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X) <- gene.names

organism <- as.character(Meta(gse)$sample_organism)

# KNN imputation
if (ncol(X) == 2) {
  X <- X[complete.cases(X), ] # KNN does not work when there are only 2 samples
} else {
  X <- X[rowSums(is.na(X)) != ncol(X), ] # remove rows with missing data
}

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

if (! is.na(argv$outrdata)) {
  save(X, pData, gene.names, organism, file = argv$outrdata)
}
