#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript overview.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.rData --rundir ~/Desktop/dgea/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --analyse "Boxplot,Volcano,PCA,Heatmap,Clustering" --distance "euclidean" --clustering "average" --dev TRUE
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

# Clustering
parser <- add_argument(parser, "--distance",
                       help = "Distance measurement methods")
parser <- add_argument(parser, "--clustering",
                       help = "HCA clustering methods")

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
factor.type     <- argv$factor
population1     <- unlist(strsplit(argv$popA, ","))
population2     <- unlist(strsplit(argv$popB, ","))
pop.name1       <- argv$popname1
pop.name2       <- argv$popname2
pop.colour1     <- "#b71c1c"  # Red
pop.colour2     <- "#0d47a1"  # Blue

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
  qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  logc <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return (logc)
}

# Boxplot
samples.boxplot <- function(data, pop.colours, pop.names, path){
  boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = Groups), outlier.shape = NA) + theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "right")+ labs(x = "Samples", y = "Expression Levels") + scale_color_manual(name = "Groups", values = pop.colours, labels = pop.names)
  # compute lower and upper whiskers
  ylim1 = boxplot.stats(data$value)$stats[c(1, 5)]

  # scale y limits based on ylim1
  boxplot <- boxplot + coord_cartesian(ylim = ylim1*1.05)

  filename <- paste(path, "boxplot.png", sep = "")
  ggsave(filename, plot = boxplot, width = 8, height = 4)

  if(isdebug){
    print("Boxplot has been produced")
  }
}

# Calculate Outliers Probabilities/ Dissimilarities
outlier.probability <- function(X, dist.method = "euclidean", clust.method = "average"){
  # Rank outliers using distance and clustering parameters
  o <- outliers.ranking(t(X),test.data = NULL, method.pars = NULL,
                        method = "sizeDiff", # Outlier finding method
                        clus = list(dist = dist.method,
                                    alg  = "hclust",
                                    meth = clust.method))
  if (isdebug) { print("Outliers have been identified") }
  return(o$prob.outliers)
}

# TODO: Add Documentation
# Clustering dendogram
clustering <- function(X, dist.method = "euclidean", clust.method = "average", exp){
  dendo  <-  hclust(dist(t(X), method = dist.method), method = clust.method)

  # Factor types
  factor <- as.factor(exp[,"factor.type"])
  population <- as.factor(exp[,"population.colour"])

  names(factor) <- exp$Sample
  names(population ) <- exp$Sample

  outliers <- outlier.probability(X, dist.method, clust.method)

  factor.cmap <- makecmap(as.numeric(factor), n = length(levels(factor)), colFn = colorRampPalette(c('black', 'green')))
  population.cmap <- makecmap(as.numeric(factor), n = length(levels(population)), colFn = colorRampPalette(c('black', 'blue')))
  outliers.cmap <- makecmap(outliers, n = 10, colFn = colorRampPalette(c('black', 'red')))

  matrix <- data.frame(Factor =  factor, Groups = population,
                       Dissimilarity = cmap( outliers, outliers.cmap))

  jColors <- with(matrix, data.frame(factors = levels(Factor),
                                     color = I(brewer.pal(nlevels(Factor), name = 'Dark2'))))

  matrix <- within(matrix,{
    Factor = jColors$color[matrix$Factor]
  })

  filename <- paste(run.dir, "clustering.png", sep = "")
  CairoPNG(file = filename, width = 1200, height = 700, xlab = "Samples")

  factor.cmap$colors         <- jColors$color
  factor.cmap$breaks         <- jColors$factors
  factor.cmap$include.lowest <- TRUE

  population.cmap$colors         <- levels(population)
  population.cmap$breaks         <- c("Group1","Group2","")
  population.cmap$include.lowest <- TRUE

  par(mar = c(6.5,6,4,3)+0.1)  # make space for color keys
  dendromat(dendo, matrix, height = 0.3, ylab = 'Distance')

  vkey(factor.cmap, 'Factors', y = 0.9, stretch = 3 )
  vkey(population.cmap, 'Groups', y = 0.6, stretch = 3)
  vkey(outliers.cmap, 'Dissimilarity', y = 0.0, stretch =2)

  dev.off()

  if (isdebug) {
    print(paste("Clustering has been performed", "with distance method:",
                argv$distance, "and clustering method:", argv$clustering))
  }
}

# Principal Component Analysis
get.pcdata <- function(Xpca){
  s <- summary(Xpca)

  exp.var <- s$importance[2, ] * 100 # Explained Variance in percentages
  cum.var <- s$importance[3, ] * 100 # Cumulative Variance in percentages
  pcnames <- names(exp.var)  # PC names

  names(exp.var) <- NULL
  names(cum.var) <- NULL

  results <- list(pcnames = pcnames, expVar = exp.var, cumVar = cum.var)

  return(results)
}

get.pcplotdata <- function(Xpca, populations){
  Xscores <- Xpca$x

  rownames(Xscores) <- NULL
  cols <- colnames(Xscores)

  # TODO: Add Documentation
  Xscores <- lapply(1:nrow(Xscores),
                    function(y) split(Xscores[, y], populations))
  names(Xscores) <- cols

  if (isdebug) { print("PCA has been calculated") }
  return(unlist(Xscores, recursive = FALSE))
}


#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################

if(isdebug){
  print("GeoDiver is starting")
  print("Libraries have been loaded")
}

if (file.exists(dbrdata)){
  load(file = dbrdata)
  if (isdebug) { print("Dataset has been loaded") }
} else {
  if (is.na(argv$geodbpath)) {
    # Load data from downloaded file
    gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE)
  } else {
    # Automatically Load GEO dataset
    gse <- getGEO(argv$accession, GSEMatrix = TRUE)
  }
  # Convert into ExpressionSet Object
  eset <- GDS2eSet(gse, do.log2 = FALSE)
}

X <- exprs(eset)  # Get Expression Data

# If not log transformed, do the log2 transformed
if (scalable(X)) {
  X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
  X <- log2(X)
}

if(isdebug){
  print(paste("Analyzing the factor", factor.type))
  print(paste("for", pop.name1,":", argv$popA))
  print(paste("against", pop.name2,":", argv$popB))
}

#############################################################################
#                        Two Population Preparation                         #
#############################################################################
# Store gene names
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names

# Phenotype Selection
pclass           <- pData(eset)[factor.type]
colnames(pclass) <- "factor.type"

# Create a data frame with the factors
expression.info  <- data.frame(pclass, Sample = rownames(pclass),
                               row.names = rownames(pclass))

# Introduce two columns to expression.info :
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
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

if (isdebug) { print("Factors and Populations have been set") }

#############################################################################
#                        Function Calling                                 #
#############################################################################

json.list <- list()

if ("Boxplot" %in% analysis.list){
  samples.boxplot(data, c(pop.colour1, pop.colour2),
                  c(pop.name1, pop.name2), path = run.dir)
}

if ("PCA" %in% analysis.list){
  Xpca <- prcomp(t(X), scale = TRUE)

  # PC individual and cumulative values
  pcdata <- get.pcdata(Xpca)
  json.list <- append(json.list, list(pc = pcdata))

  # PC scatter plot
  pcplotdata <- get.pcplotdata(Xpca, expression.info[, "population"])

  # adding both data to the json list
  json.list <- append(json.list, list(pcdata = pcplotdata))
}

if ("Clustering" %in% analysis.list){
  clustering(X, dist.method, clust.method, expression.info)
}

if (length(json.list) != 0) {
  filename <- paste(run.dir, "data.json", sep = "")
  write(toJSON(json.list, digits=I(4)), filename )
}
