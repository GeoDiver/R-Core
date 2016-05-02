#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript dgea.R --dbrdata ~/Desktop/GDS5093.rData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --analyse "Volcano,PCA,Heatmap" --topgenecount 250 --foldchange 0.0 --thresholdvalue 0.005 --distance "euclidean" --clustering "average" --clusterby "Complete" --heatmaprows 100 --dendrow --dendcol --adjmethod fdr --dev
# ---------------------------------------------------------#

#############################################################################
#                        Import Libraries                                   #
#############################################################################

# load required libraries and silence library loading messages on command line
suppressMessages(library("argparser"))     # Argument passing
suppressMessages(library("Cairo"))         # Plots saving
suppressMessages(library("dendextend"))    # Dendogram extended functionalities
suppressMessages(library("DMwR"))          # Outlier Prediction for clustering
suppressMessages(library("GEOquery"))      # GEO dataset Retrieval
suppressMessages(library("ggplot2"))       # Graphs designing
suppressMessages(library("jsonlite"))      # Convert R object to JSON format
suppressMessages(library("limma"))         # Differencial Gene Expression Analysis
suppressMessages(library("pheatmap"))      # Heatmap Generating
suppressMessages(library("plyr"))          # Splitting, Applying and Combining Data
suppressMessages(library("RColorBrewer"))  # Import Colour Pallete
suppressMessages(library("reshape2"))      # Prepare dataset for ggplot
suppressMessages(library("squash"))        # Clustering Dendogram

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--rundir",
                       help = "The outout directory where graphs get saved")
parser <- add_argument(parser, "--dbrdata",
                       help = "Downloaded GEO dataset full path")
parser <- add_argument(parser, "--analyse", nargs = "+",
                       help = "List of analysis to be performed")
parser <- add_argument(parser, "--geodbpath",
                       help = "GEO Dataset full path")
parser <- add_argument(parser, "--dev", flag = TRUE, short = '-d',
                       help = "verbose version")

# Sample Parameters
parser <- add_argument(parser, "--factor",
                       help = "Factor type to be classified by")
parser <- add_argument(parser, "--popA", nargs = "+",
                       help = "Group A - all the selected phenotypes (at least one)")
parser <- add_argument(parser, "--popB", nargs = "+",
                       help = "Group B - all the selected phenotypes (at least one)")
parser <- add_argument(parser, "--popname1",
                       help = "name for Group A")
parser <- add_argument(parser, "--popname2",
                       help = "name for Group B")

# Toptable
parser <- add_argument(parser, "--topgenecount",
                       help = "Number of top genes to be used")
parser <- add_argument(parser, "--adjmethod",
                       help = "Adjust P-values for Multiple Comparisons")

# Heatmap
parser <- add_argument(parser, "--heatmaprows",
                       help = "Number of genes show in the heatmap")
parser <- add_argument(parser, "--dendrow", flag = TRUE,
                       help = "Flag to display dendogram for Genes")
parser <- add_argument(parser, "--dendcol", flag = TRUE,
                       help = "Flag to display dendogram for Samples")
parser <- add_argument(parser, "--clusterby",
                       help = "Cluster by complete dataset or toptable")
# Clustering
parser <- add_argument(parser, "--distance",
                       help = "Distance measurement methods")
parser <- add_argument(parser, "--clustering",
                       help = "HCA clustering methods")

# Volcano plot Parameters
parser <- add_argument(parser, "--foldchange",
                       help = "fold change cut off")
parser <- add_argument(parser, "--thresholdvalue",
                       help = "threshold value cut off")

# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                        Command Line Arguments Retrieval                   #
#############################################################################
split_arg <- function(vector_arg) {
  pos <- unlist(gregexpr("[^\\\\],", vector_arg, perl=TRUE))
  vect <- substring(vector_arg, c(1, pos+2), c(pos, nchar(vector_arg)))
  vect <- gsub("\\\\,", ",", vect) # replace \\, with ,
  vect <- gsub("\\\"", '"', vect)  # replace \" with "
  vect <- gsub('\\\\-', '-', vect) # replace \\- with -
  return (vect)
}

# General Parameters
run.dir       <- argv$rundir
dbrdata       <- argv$dbrdata
analysis.list <- split_arg(argv$analyse)
isdebug       <- argv$dev

# Sample Parameters
factor.type   <- argv$factor
population1   <- split_arg(argv$popA)
population2   <- split_arg(argv$popB)
pop.name1     <- argv$popname1
pop.name2     <- argv$popname2
pop.colour1   <- "#b71c1c"  # Red
pop.colour2   <- "#0d47a1"  # Blue

# Toptable
topgene.count   <- as.numeric(argv$topgenecount)
toptable.sortby <- "p"
adjmethod_opt   <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                     "fdr", "none")
adj.method <- ifelse(argv$adjmethod %in% adjmethod_opt, argv$adjmethod, "fdr")

# Heatmap
heatmap.rows <- as.numeric(argv$heatmaprows)
dendrow      <- as.logical(argv$dendrow)
dendcol      <- as.logical(argv$dendcol)
cluster.by   <- argv$clusterby

# Clustering
dist.opt     <- c("euclidean", "maximum", "manhattan", "canberra", "binary",
                  "minkowski")
dist.method  <- ifelse(argv$distance %in% dist.opt, argv$distance, "euclidean")

clust.opt    <- c("ward.D", "ward.D2", "single", "complete", "average",
                  "mcquitty", "median", "centroid")
clust.method <- ifelse(argv$clustering %in% clust.opt, argv$clustering,
                       "average")

# Volcano plot Parameters
fold.change     <- as.numeric(argv$foldchange)
threshold.value <- as.numeric(argv$thresholdvalue)

#############################################################################
#                          Load Functions                                   #
#############################################################################

# Check if the run directory exists and if not, create directory...
check.run.dir <- function(run.dir) {
  if (!dir.exists(file.path(run.dir))) {
    dir.create(file.path(run.dir))
  }
}

# Calculate Outliers Probabilities/ Dissimilarities
outlier.probability <- function(X, dist.method = "euclidean",
                                clust.method = "average") {
  # Rank outliers using distance and clustering parameters
  o <- outliers.ranking(t(X),test.data = NULL, method.pars = NULL,
                        method = "sizeDiff", # Outlier finding method
                        clus = list(dist = dist.method,
                                    alg  = "hclust",
                                    meth = clust.method))
  if (isdebug) print("DGEA: Outliers have been identified")
  return(o$prob.outliers)
}

find.toptable <- function(X, newpclass, toptable.sortby, topgene.count) {
  design  <- model.matrix(~0 + newpclass)  # creates a design (or model) matrix

  # plots linear model for each gene and estimate fold changes + standard errors
  fit     <- lmFit(X, design)

  # set contrasts for two groups
  contrasts <- makeContrasts(contrasts = "newpclassGroup1-newpclassGroup2",
                             levels = design)

  fit <- contrasts.fit(fit, contrasts)

  # empirical bayes smoothing to standard errors
  fit <- eBayes(fit, proportion = 0.01)

  # create top Table
  toptable <- topTable(fit, adjust.method = adj.method,
                       sort.by = toptable.sortby, number = topgene.count)
  if (isdebug) {
    print(paste("DGEA: TopTable has been produced",
          "for", topgene.count, "genes with the cut-off method:", adj.method))
  }
  return(toptable)
}

# Heatmap
heatmap <- function(X.matrix, X, exp, heatmap.rows = 100, dendogram.row, dendogram.col,
                    dist.method, clust.method, path) {
  col.pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlGn")))(100)

  # Annotation column for samples
  ann.col <- data.frame(Population = exp[, "population"])

   # If there are only few factor levels, then only show colour annotations
  if (length(levels(exp[, "factor.type"])) < 10) {
    ann.col$Factor    <- exp[, "factor.type"]
    colnames(ann.col) <- c("Population", factor.type)
  }

  # Clustering based on complete dataset or only toptable data
  if (cluster.by == "Complete") {
    exp.data <- X
  } else {
    exp.data <- X.matrix
  }

  # Column dendogram
  if (dendogram.col != TRUE) {
    hc <- FALSE
    # Keep a gap between two groups
    column.gap <- length( (which(ann.col[, "Population"] == "Group1") == T) )
  } else {
    # calculate heirachical clustering
    hc <- hclust(dist(t(exp.data), method = dist.method), method = clust.method)

    # Find outlier ranking/ probability
    outliers <- outlier.probability(exp.data, dist.method, clust.method)

    # Add dissimilarity annotation to the samples annotation
    ann.col$Dissimilarity <- outliers
    column.gap <- 0
  }

  rownames(ann.col) <- exp[, "Sample"]

  # Limit no of heatmap rows
  if (nrow(X.matrix) < heatmap.rows) {
    hdata <- X.matrix
  } else {
    hdata <- X.matrix[1:heatmap.rows, ]
  }

  filename <- file.path(path, "dgea_heatmap.svg")
  CairoSVG(file = filename)

  pheatmap(hdata,
           cluster_row    = dendogram.row,
           cluster_cols   = hc,
           annotation_col = ann.col,
           legend         = TRUE,
           color          = col.pal,
           fontsize       = 6.5,
           fontsize_row   = 3.0,
           fontsize_col   = 3.5,
           gaps_col       = column.gap)
  dev.off()

  if (isdebug) {
    print(paste("DGEA: Heatmap has been created"))
    if (dendrow==TRUE) print("DGEA: with a dendogram for rows")
    if (dendcol==TRUE) print("DGEA: and a dendogram for columns")
  }
}

# Apply Bonferroni cut-off as the default thresold value
# fold.change and threshold value are no longer used to find significant gene. 
# But the necessary code is left due to plans to extend the functionality with
# those thwo parameters.
volcanoplot <- function(toptable, fold.change, t = 0.05 / length(gene.names),
                        path) {
  # Select only genes which are in toptable
  toptable$Significant <- c(rep(TRUE,topgene.count),
                            rep(FALSE,length(toptable$ID) - topgene.count))

  # Construct the plot object
  vol <- ggplot(data = toptable, aes(x = toptable$logFC, y = -log10(toptable$P.Value), colour = Significant)) +
      geom_point(alpha = 0.4, size = 1.75)  + xlim(c(-max(toptable$logFC) - 0.1, max(toptable$logFC) + 0.1)) + ylim(c(0, max(-log10(toptable$P.Value)) + 0.5)) + xlab("log2 fold change") + ylab("-log10 p-value")

  # File saving as png
  filename <- file.path(path, "dgea_volcano.png")
  ggsave(filename, plot = vol, height = 6, width = 6)

  if (isdebug) {
    print(paste("DGEA: Volcanoplot has been produced",
          "with a foldchange of:", fold.change, "and threshold of:", threshold.value))
  }
}

get.volcanodata <- function(toptable) {
  vol.list <- list(genes = toptable$ID,
                   logFC = round(toptable$logFC, 3),
                   pVal  = -log10(toptable$P.Value))
  return(vol.list)
}

#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################

if (file.exists(dbrdata)) {
  if (isdebug) print("DGEA: Loading Database data.")
  load(file = dbrdata)
} else {
  cat("ERROR: Data input error. Provide valid GDS dataset!", file=stderr())
  quit(save = "no", status = 200, runLast = FALSE)
}

check.run.dir(run.dir)

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Phenotype Selection
pclass           <- pData[factor.type]
colnames(pclass) <- "factor.type"

# Create a data frame with the factors
expression.info  <- data.frame(pclass, Sample = rownames(pclass),
                               row.names = rownames(pclass))

# Introduce two columns to expression.info :
#   1. population - new two groups/populations, NA for unselected samples
#   2. population.colour - colour for two new two groups, black for unselected samples
expression.info <- within(expression.info, {
  population        <- ifelse(factor.type %in% population1, "Group1",
                         ifelse(factor.type %in% population2, "Group2", NA))
  population.colour <- ifelse(factor.type %in% population1, pop.colour1,
                        ifelse(factor.type %in% population2, pop.colour2,
                        "#000000"))
})

# Convert population column to a factor
expression.info$population <- as.factor(expression.info$population)

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info), ]
X <- X[, (colnames(X) %in% rownames(expression.info))]

# Data preparation for ggplot-Boxplot
data <- within(melt(X), {
  phenotypes <- expression.info[Var2, "factor.type"]
  Groups     <- expression.info[Var2, "population.colour"]
})

# Created a new Phenotype class
newpclass        <- expression.info$population
names(newpclass) <- expression.info$Sample

if (isdebug) print("DGEA: Factors and Populations have been set")

#############################################################################
#                        Function Calling                                   #
#############################################################################

# empty list to collect all data need to be displayed in plotly
json.list <- list()

# Toptable
toptable <- find.toptable(X, newpclass, toptable.sortby, topgene.count)

# Adding to JSON file
temp.toptable <- toptable
names(temp.toptable) <- NULL
json.list <- append(json.list, list(tops = temp.toptable))

# Filter toptable data from X
X.toptable <- X[as.numeric(rownames(toptable)), ]

# save toptable expression data
filename <- file.path(run.dir,"dgea_toptable.RData")
save(X.toptable, expression.info, file = filename)

# save tab delimited
filename <- file.path(run.dir, "dgea_toptable.tsv")
write.table(toptable, filename, col.names=NA, sep = "\t" )

if (isdebug) print(paste("DGEA: Analysis to be performed:", argv$analyse))

if ("Volcano" %in% analysis.list) {
  # Get data for volcanoplot
  toptable.all <- find.toptable(X, newpclass, toptable.sortby,
                                length(gene.names))
  # Draw volcano plit
  volcanoplot(toptable.all, fold.change, threshold.value, run.dir)

  # save volcanoplot top data as JSON
  volcanoplot.data <- get.volcanodata(toptable)
  json.list        <- append(json.list, list(vol = volcanoplot.data))
}

if ("Heatmap" %in% analysis.list) {
  heatmap(X.toptable,X, expression.info, heatmap.rows = heatmap.rows, dendrow,
          dendcol, dist.method, clust.method, run.dir)
}

if (length(json.list) != 0) {
  filename <- file.path(run.dir, "dgea_data.json")
  write(toJSON(json.list, digits=I(4)), filename) # 4 decimal places
}
