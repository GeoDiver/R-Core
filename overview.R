#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript overview.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.rData --rundir ~/Desktop/dgea/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --analyse "Boxplot,PCA" --dev TRUE
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
pop.colour1     <- "#e199ff" # Purple   
pop.colour2     <- "#96ca00" # Green

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
  #  produce sample quantiles corresponding to the given probabilities
  qx <- as.numeric(quantile(X, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  logc <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return (logc)
}

# Boxplot
samples.boxplot <- function(data, pop.colours, pop.names, path){
 
  boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = Groups), outlier.shape = NA) + theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "right")+ labs(x = "Samples", y = "Expression Levels") + scale_color_manual(name = "Groups", values = pop.colours, labels = pop.names)
  
  # compute lower and upper whiskers to set the y axis margin
  ylim1 = boxplot.stats(data$value)$stats[c(1, 5)]

  # scale y limits based on ylim1
  boxplot <- boxplot + coord_cartesian(ylim = ylim1*1.05)

  filename <- paste(path, "boxplot.png", sep = "")
  ggsave(filename, plot = boxplot, width = 8, height = 4)

  if (isdebug){
    print("Overview: Boxplot has been produced")
  }
}

# Principal Component Analysis
get.pcdata <- function(Xpca){
  s <- summary(Xpca)

  exp.var <- s$importance[2, ] * 100 # Explained Variance in percentages
  cum.var <- s$importance[3, ] * 100 # Cumulative Variance in percentages
  pcnames <- names(exp.var)          # PC names

  names(exp.var) <- NULL
  names(cum.var) <- NULL

  results <- list(pcnames = pcnames, expVar = exp.var, cumVar = cum.var)

  return(results)
}

get.pcplotdata <- function(Xpca, populations){
  
  Xscores <- Xpca$x

  sample.names <- rownames(Xscores)
  rownames(Xscores) <- NULL
  cols <- colnames(Xscores)

  # Take each column of Xscore (as temp variable y) and split based on the population
  Xscores <- lapply(1:nrow(Xscores),
                    function(y) split(Xscores[, y], populations))
  names(Xscores) <- cols
  
  # Unlist them but not to the dept. outer most list has unlisted
  pc <- unlist(Xscores, recursive = FALSE)
  
  # Split sample names by population and add to the final list
  complete.data <- c(split(sample.names, populations),pc)
    
  if (isdebug) { print("Overview: PCA has been calculated") }
  return(complete.data)
}


#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################

if(isdebug){
  print("Overview: GeoDiver is starting")
  print("Overview: Libraries have been loaded")
}

if (file.exists(dbrdata)){
  load(file = dbrdata)
} else {
  if (is.na(argv$geodbpath)) {
    # Load data from downloaded file
    gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE)
  } else {
    # Automatically Load GEO dataset
    gse <- getGEO(accession, GSEMatrix = TRUE)
  }
  # Convert into ExpressionSet Object
  eset <- GDS2eSet(gse, do.log2 = FALSE)
}

#############################################################################
#                           Data Preprocessing                              #
#############################################################################

X <- exprs(eset)  # Get Expression Data

# Remove NA Data from the dataset
not.null.indexes <- which(complete.cases(X[,])==TRUE)
X <- X[not.null.indexes,]

# If not log transformed, do the log2 transformed
if (scalable(X)) {
  X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
  X <- log2(X)
}

if (isdebug){print("Overview: Data Preprocessed!")}

#############################################################################
#                        Two Population Preparation                         #
#############################################################################
# Store gene names
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names[not.null.indexes]

# Phenotype selection
pclass           <- pData(eset)[factor.type]
colnames(pclass) <- "factor.type"

# Create a data frame with the factors
expression.info  <- data.frame(pclass, Sample = rownames(pclass),
                               row.names = rownames(pclass))

# Introduce two columns to expression.info :
#   1. population - new two groups, if not NA
#   2. population.colour - colour for two new two groups, if not black colour
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

if (isdebug) { print("Overview: Factors and Populations have been set") }

#############################################################################
#                        Function Calling                                 #
#############################################################################

json.list <- list()

if ("Boxplot" %in% analysis.list){
  samples.boxplot(data, c(pop.colour2, pop.colour1),
                  c(pop.name2, pop.name1), path = run.dir)
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

if (length(json.list) != 0) {
  # Write to a json file with 4 decimal places
  filename <- paste(run.dir, "data.json", sep = "")
  write(toJSON(json.list, digits=I(4)), filename )
}
