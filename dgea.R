#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript dgea.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.rData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --analyse "Volcano,PCA,Heatmap" --topgenecount 250 --foldchange 0.0 --thresholdvalue 0.005 --distance "euclidean" --clustering "average" --clusterby "Complete" --heatmaprows 100 --dendrow TRUE --dendcol TRUE --adjmethod fdr --dev TRUE
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
parser <- add_argument(parser, "--dev",
                       help = "The output directory where graphs get saved")

# Sample Parameters
parser <- add_argument(parser, "--accession",
                       help = "Accession Number of the GEO Database")
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
parser <- add_argument(parser, "--dendrow",
                       help = "Boolean value for display dendogram for Genes")
parser <- add_argument(parser, "--dendcol",
                       help = "Boolean value for display dendogram for Samples")
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

# General Parameters
run.dir         <- argv$rundir
dbrdata         <- argv$dbrdata
analysis.list   <- unlist(strsplit(argv$analyse, ","))

# Sample Parameters
accession       <- argv$accession
factor.type     <- argv$factor
population1     <- unlist(strsplit(argv$popA, ","))
population2     <- unlist(strsplit(argv$popB, ","))
pop.name1       <- argv$popname1
pop.name2       <- argv$popname2
pop.colour1     <- "#b71c1c"  # Red
pop.colour2     <- "#0d47a1"  # Blue

# Toptable
topgene.count     <- as.numeric(argv$topgenecount)
toptable.sortby   <- "p"
adjmethod_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                       "fdr", "none")
if (argv$adjmethod %in% adjmethod_options) {
  adj.method <- argv$adjmethod
} else {
  adj.method <- "fdr"
}

# Heatmap
heatmap.rows <- as.numeric(argv$heatmaprows)
dendrow      <- as.logical(argv$dendrow)
dendcol      <- as.logical(argv$dendcol)
cluster.by   <- argv$clusterby

# Clustering
distance_options <- c("euclidean", "maximum", "manhattan", "canberra",
                      "binary", "minkowski")
if (argv$distance %in% distance_options){
  dist.method <- argv$distance
} else {
  dist.method <- "euclidean"
}

clustering_options <- c("ward.D", "ward.D2", "single", "complete", "average",
                        "mcquitty", "median", "centroid")
if (argv$clustering %in% clustering_options){
  clust.method <- argv$clustering
} else {
  clust.method <- "average"
}

# Volcano plot Parameters
fold.change     <- as.numeric(argv$foldchange)
threshold.value <- as.numeric(argv$thresholdvalue)

if (!is.na(argv$dev)) {
  isdebug <- argv$dev
} else {
  isdebug <- FALSE
}

#############################################################################
#                          Load Functions                                   #
#############################################################################

# auto-detect if data is log transformed
scalable <- function(X) {
  # sample quantiles corresponding to the given probabilities
  qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  logc <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return (logc)
}

# Calculate Outliers Probabilities/ Dissimilarities
outlier.probability <- function(X, dist.method = "euclidean", clust.method = "average"){

 # Rank outliers using distance and clustering parameters
  o <- outliers.ranking(t(X),test.data = NULL, method.pars = NULL,
                        method = "sizeDiff", # Outlier finding method
                        clus = list(dist = dist.method,
                                    alg  = "hclust",
                                    meth = clust.method))
  if (isdebug) print("DGEA: Outliers have been identified")
  return(o$prob.outliers)
}

find.toptable <- function(X, newpclass, toptable.sortby, topgene.count){

  # creates a design (or model) matrix
  design  <- model.matrix(~0 + newpclass)

  # plots linear model for each gene and estimate fold changes and standard errors
  fit     <- lmFit(X, design)

  # set contrasts for two groups
  contrasts <- makeContrasts(contrasts = "newpclassGroup1-newpclassGroup2",
                             levels = design)

  fit <- contrasts.fit(fit, contrasts)

  # empirical bayes smoothing to standard errors
  fit <- eBayes(fit, proportion = 0.01)

  # create top Table
  toptable <- topTable(fit, adjust.method = adj.method, sort.by = toptable.sortby,
                       number = topgene.count)
  if(isdebug){
    print(paste("DGEA: TopTable has been produced",
          "for", topgene.count, "genes with the cut-off method:", adj.method))
  }
  return(toptable)
}

# Heatmap
heatmap <- function(X.matix, X, exp, heatmap.rows = 100, dendogram.row, dendogram.col,
                    dist.method, clust.method, path){

  col.pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlGn")))(100)

  # Annotation column for samples
  ann.col <- data.frame(Population = exp[, "population"])

   # If there are only few factor levels, then only show colour annotations
  if(length(levels(exp[, "factor.type"])) < 10){
      ann.col$Factor <- exp[, "factor.type"]
      colnames(ann.col) <- c("Population", factor.type)
  }

  # Clustering based on complete dataset or only toptable data
  if (cluster.by == "Complete"){
      exp.data <- X
  }else{
      exp.data <- X.matix
  }

  # Column dendogram
  if (dendogram.col == TRUE){

    # calculate heirachical clustering
    hc <- hclust(dist(t(exp.data), method = dist.method), method = clust.method)

    # Find outlier ranking/ probability
    outliers <- outlier.probability(exp.data, dist.method, clust.method)

    # Add dissimilarity annotation to the samples annotation
    ann.col$Dissimilarity <- outliers

    column.gap <- 0

  } else {

    hc <- FALSE

    # Keep a gap between two groups
    column.gap <- length( (which(ann.col[, "Population"] == "Group1") == T) )
  }

  rownames(ann.col) <- exp[, "Sample"]

  # Limit no of heatmap rows
  if(nrow(X.matix) < heatmap.rows){
      hdata <- X.matix                    # Show all
  }else{
      hdata <- X.matix[1:heatmap.rows, ]  # Limit to user specified limit
  }

  filename <- paste(path, "dgea_heatmap.svg", sep = "")
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
# fold.change and threshold value are not used to find significant gene. But it kept
# remain as there is a plan to extend the functionality with those two parameters.
volcanoplot <- function(toptable, fold.change, t = 0.05 / length(gene.names), path){

  # Select only genes which are in toptable
  toptable$Significant <- c(rep(TRUE,topgene.count),
                            rep(FALSE,length(toptable$ID) - topgene.count))

  # Construct the plot object
  vol <- ggplot(data = toptable, aes(x = toptable$logFC, y = -log10(toptable$P.Value), colour = Significant)) +
      geom_point(alpha = 0.4, size = 1.75)  + xlim(c(-max(toptable$logFC) - 0.1, max(toptable$logFC) + 0.1)) + ylim(c(0, max(-log10(toptable$P.Value)) + 0.5)) + xlab("log2 fold change") + ylab("-log10 p-value")

  # File saving as png
  filename <- paste(path, "dgea_volcano.png", sep = "")
  ggsave(filename, plot = vol, height = 6, width = 6)

  if(isdebug){
    print(paste("DGEA: Volcanoplot has been produced",
          "for a foldchange of:", fold.change, "and threshold of:", threshold.value))
  }
}

get.volcanodata <- function(toptable){
  vol.list <- list( genes = toptable$ID,
                    logFC = round(toptable$logFC, 3),
                    pVal  = -log10(toptable$P.Value))
  return(vol.list)
}

#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################

if (isdebug){
  print("DGEA: GeoDiver is starting")
  print("DGEA: Libraries have been loaded")
}

if (file.exists(dbrdata)){
  load(file = dbrdata)
  if (isdebug) print("DGEA: Dataset has been loaded")
} else {
  tryCatch({
    # Automatically Load GEO dataset
    gse <- getGEO(accession, GSEMatrix = TRUE)

    # Convert into ExpressionSet Object
    eset <- GDS2eSet(gse, do.log2 = FALSE)
  },error=function(e){
      print("ERROR:Data input error. Provide valid GDS dataset!")
      # Exit with error code 1
      quit(save = "no", status = 1, runLast = FALSE)
  })
}

#############################################################################
#                           Data Preprocessing                              #
#############################################################################

X <- exprs(eset)  # Get Expression Data

# Remove NA data from the dataset
not.null.indexes <- which(complete.cases(X[,])==TRUE)
X <- X[not.null.indexes,]

# If not log transformed, do the log2 transformed
if (scalable(X)) {
    X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
    X <- log2(X)
}

if (isdebug) print("DGEA: Data Preprocessed!")

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Store gene names
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names[not.null.indexes]

# Phenotype Selection
pclass           <- pData(eset)[factor.type]
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
newpclass           <- expression.info$population
names(newpclass)    <- expression.info$Sample

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
filename <- paste(run.dir,"dgea_toptable.RData", sep = "")
save(X.toptable, expression.info, file = filename)

# save tab delimited
filename <- paste(run.dir, "dgea_toptable.tsv", sep = "")
write.table(toptable, filename, col.names=NA, sep = "\t" )

if(isdebug) print(paste("DGEA: Analysis to be performed:", argv$analyse))


if ("Volcano" %in% analysis.list){
    # Get data for volcanoplot
    toptable.all <- find.toptable(X, newpclass, toptable.sortby,
                                  length(gene.names))
    # Draw volcano plit
    volcanoplot(toptable.all, fold.change, threshold.value, run.dir)

    # save volcanoplot top data as JSON
    volcanoplot.data <- get.volcanodata(toptable)
    json.list <- append(json.list, list(vol = volcanoplot.data))
}

if ("Heatmap" %in% analysis.list){

    heatmap(X.toptable,X, expression.info, heatmap.rows = heatmap.rows,
            dendrow, dendcol, dist.method, clust.method, run.dir)
}

if (length(json.list) != 0){
    # Write to a json file with 4 decimal places
    filename <- paste(run.dir, "dgea_data.json", sep = "")
    write(toJSON(json.list, digits=I(4)), filename )
}
