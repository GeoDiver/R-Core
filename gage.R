#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Gene Set Enrichment Analysis             #
# Rscript gage.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.RData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control"  --comparisontype ExpVsCtrl --genesettype KEGG --distance "euclidean" --clustering "average" --clusterby "Complete" --heatmaprows 100 --dendrow TRUE --dendcol TRUE --dev TRUE
# ---------------------------------------------------------#

#############################################################################
#                        Import Libraries                                   #
#############################################################################

# load required libraries and silence library loading messages on command line
suppressMessages(library("argparser"))     # Argument passing
suppressMessages(library("Cairo"))         # Plots saving
suppressMessages(library("DMwR"))          # Outlier Prediction for clustering
suppressMessages(library("gage"))          # Does the analysis
suppressMessages(library("gageData"))      # Lets data be used by gage
suppressMessages(library("GEOquery"))      # GEO dataset Retrieval
suppressMessages(library("GO.db"))         # Loads GO database
suppressMessages(library("jsonlite"))      # Convert R object to JSON format
suppressMessages(library("pathview"))      # Interaction networks & used to get ENTREZ IDs
suppressMessages(library("pheatmap"))      # Used to create heatmap
suppressMessages(library("RColorBrewer"))  # Color palette for heatmap
suppressMessages(library("org.Mm.eg.db"))  # Species database

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Paeameters
parser <- add_argument(parser, "--accession",
                       help = "Accession Number of the GEO Database")
parser <- add_argument(parser, "--dbrdata",
                       help = "Downloaded GEO dataset full path")
parser <- add_argument(parser, "--rundir",
                       help = "The output directory where graphs get saved")
parser <- add_argument(parser, "--dev",
                       help = "The output directory where graphs get saved")

# Sample Parameters
parser <- add_argument(parser, "--popA",
                       help = "GroupA - all the selected phenotypes (atleast one)",
                       nargs = "+")
parser <- add_argument(parser, "--popB",
                       help = "GroupB - all the selected phenotypes (atleast one)",
                       nargs = "+")
parser <- add_argument(parser, "--factor",
                       help = "Factor type to be classified by")

# Gage Parameters
parser <- add_argument(parser, "--comparisontype", help = "ExpVsCtrl or ExpVsExp")
parser <- add_argument(parser, "--genesettype",
                       help = "KEGG - KEGG Database or
                       BP - GO Biological Process or
                       MF - GO molecular function
                       or CC - Cellular Component")

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


# allows arguments to be run via the command line
argv <- parse_args(parser)

# #############################################################################
#                         Command Line Arguments Retrieval                   #
# #############################################################################

# General Parameters
rundir          <- argv$rundir
dbrdata         <- argv$dbrdata
isdebug         <- ifelse(!is.na(argv$dev), argv$dev, FALSE)

# Sample Parameters
accession   <- argv$accession
factor.type <- as.character(argv$factor)
population1 <- unlist(strsplit(argv$popA, ","))
population2 <- unlist(strsplit(argv$popB, ","))

pop.colour1     <- "#b71c1c"  # Red
pop.colour2     <- "#0d47a1"  # Blue

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

# Gage parameters
comparison.type <- argv$comparisontype  # "ExpVsCtrl" or "ExpVsExp"
geneset.type    <- argv$genesettype     # "KEGG"  or "BP" or "MF" or "CC"

#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################
if (isdebug){
    print("GAGE: GeoDiver is starting")
}

if (isdebug){
    print("GAGE: Libraries have been loaded")
}

if (file.exists(dbrdata)){
    load(file = dbrdata)
} else {
    # Automatically Load GEO dataset
    gse <- getGEO(accession, GSEMatrix = TRUE)
    
    # Convert into ExpressionSet Object
    eset <- GDS2eSet(gse, do.log2 = TRUE)
}

# Auto-detect if data is log transformed
scalable <- function(X) {
    qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    logc <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    return (logc)
}
#############################################################################
#                           Data Preprocessing                              #
#############################################################################

# Get Expression Data
X <- exprs(eset)

# Remove NA Data from the dataset
not.null.indexes <- which(complete.cases(X[, ]) == TRUE)
X <- X[not.null.indexes, ]

# If not log transformed, do the log2 transformed
if (scalable(X)) {
    X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
    X <- log2(X)
}

if (isdebug){print("GAGE: Data Preprocessed!")}

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

if (isdebug){print(paste("GAGE: Factor :", factor.type))}
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names[not.null.indexes]

# Phenotype Selection
pclass           <- pData(eset)[factor.type]
colnames(pclass) <- "factor.type"

# Create a data frame with the factors
expression.info  <- data.frame(pclass,
                               Sample = rownames(pclass),
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

# Get sample indexes and sample names
Group1 <-  which(expression.info[, "population"] == "Group1")
Group1names <- expression.info[Group1, "Sample"]

Group2 <-  which(expression.info[, "population"] == "Group2")
Group2names <- expression.info[Group2, "Sample"]

if (isdebug ){
    print("GAGE: Factors and Populations have been set")
}

#############################################################################
#                            Data Preparation                               #
#############################################################################

# Loading data on kegg packages and species
data(bods)
data(korg)

# Retrieve orgamism from metadata
organism    <- as.character(Meta(gse)$sample_organism)

# Retrieve KEGG code and package for the organism
keggcode.organism <-as.character(korg[which(korg[, "scientific.name"] == organism), "kegg.code"])
package <-as.character(bods[which(bods[, "kegg code"] == keggcode.organism), "package"])

# Create two column table containing entrez IDs for geodataset
id.map.refseq <- id2eg(ids =  gene.names, category = "SYMBOL", 
                       pkg.name = package, org = as.character(keggcode.organism))

# Replace gene symbols with ENTREZ ID in dataset matrix
rownames(X) <- id.map.refseq[, 2]

# Remove rows without ENTREZ IDs
X <- X[which(is.na(rownames(X)) == FALSE), ]

geo.dataset <- X

if(isdebug){
    print("GAGE: Data Preparation completed")
}

#############################################################################
#                          Gage  Data Loading                               #
#############################################################################

if (geneset.type == "KEGG") { # KEGG datasets
    data(kegg.gs)
    kg.org  <- kegg.gsets(organism)               # picks out orgamism gene sets
    dbdata <- kg.org$kg.sets[kg.org$sigmet.idx]
} else { # GO Datasets
    common.name <- as.character(bods[which(bods[, "kegg code"] == keggcode.organism ), "species"])
    go.hs <- go.gsets(species = common.name)      # use species column of bods
    if(geneset.type == "BP"){                     # BP = Biological Process
        dbdata <- go.hs$go.sets[go.hs$go.subs$BP] 
    } else if(geneset.type == "MF"){              # MF = molecular function
        dbdata <- go.hs$go.sets[go.hs$go.subs$MF]
    } else if(geneset.type == "CC"){              # CC = cellular component
        dbdata <- go.hs$go.sets[go.hs$go.subs$CC]
    }
}

#############################################################################
#                            Heatmap                                        #
#############################################################################

# Calculate Outliers Probabilities/ Dissimilarities
outlier.probability <- function(X, dist.method = "euclidean", clust.method = "average"){
    
    # Rank outliers using distance and clustering parameters
    o <- outliers.ranking(t(X),test.data = NULL, method.pars = NULL,
                          method = "sizeDiff", # Outlier finding method
                          clus = list(dist = dist.method,
                                      alg  = "hclust",
                                      meth = clust.method))
    if (isdebug) { print("GAGE: Outliers have been identified") }
    return(o$prob.outliers)
}

get.heatmap <- function(analysis.stats, analysis.type){
    
    analysis.heatmap <- t(analysis.stats)
    analysis.heatmap <- analysis.heatmap
    row.names(analysis.heatmap) <- gsub("(stats.)", "", row.names(analysis.heatmap))
    
    col.pal <- colorRampPalette(rev(
        RColorBrewer::brewer.pal(11, "RdYlGn")))(100)
    
    if (analysis.type == "ExpVsCtrl"){ # Remove control group (Group 2)
        X.mat <- X[,Group1]
        exp.info <- expression.info[expression.info[,"population"]== "Group1", ]
    } else {
        X.mat <- X
        exp.info <- expression.info
    }
    
    # Column dendogram
    if (dendcol == TRUE){
        # calculate heirachical clustering
        hc <- hclust(dist(t(X.mat), method = dist.method), method = clust.method)
        
        # Find outlier ranking/ probability
        outliers <- outlier.probability(X.mat, dist.method, clust.method)
        
        # Annotation columns
        ann.col <- data.frame(Population    = exp.info[, "population"],
                              Factor        = exp.info[, "factor.type"],
                              Dissimilarity = outliers)
        colnames(ann.col) <- c("Population", factor.type,"Dissimilarity")
    } else {
        hc <- FALSE
        # Annotation columns
        ann.col <- data.frame(Population    = exp.info[, "population"],
                              Factor        = exp.info[, "factor.type"])
        colnames(ann.col) <- c("Population", factor.type)
    }
    trans.analysis <- t(analysis.heatmap)
    if (nrow(trans.analysis) < heatmap.rows){
        hdata <- trans.analysis                    # Show all 
    } else {
        hdata <- trans.analysis[1:heatmap.rows, ]  # Limit to user specified limit
    }
    
    filename <- paste(rundir, "gage_heatmap.svg", sep="")
    CairoSVG(file = filename)
    pheatmap(hdata,
             cluster_row    = dendrow,
             cluster_cols   = hc,
             annotation_col = ann.col,
             color          = col.pal,
             fontsize       = 6.5,
             fontsize_row   = 3.0,
             fontsize_col   = 3.5)
    dev.off()
    
    if(isdebug ){
        print(paste("GAGE: Saved heatmap", filename))
    }
}

#############################################################################
#                         GAGE ANALYSIS                                     #
#############################################################################

gage.analysis <- function(set.type, analysis.type = "ExpVsCtrl", ref.group = G2, samp.group = G1, compare.option = "unpaired"){
    
    analysis <- gage(geo.dataset, gsets = set.type,
                     ref = G2, samp = G1,
                     same.dir = F, compare = compare.option )
    
    # Returns number of two-direction significantly enriched gene sets
    analysis.sig <- sigGeneSet(analysis)
    
    if (nrow(analysis.sig$greater) > 0){
        # Formatting and preparation for heatmap
        analysis.sig <- as.data.frame(analysis.sig)
        analysis.stats <- analysis.sig[, grep("^stats.GSM", names(analysis.sig), value = TRUE)]
        
        # Split each pathway names into pathway ID and pathway name
        m <- regmatches(rownames(analysis.sig), regexpr(" ", rownames(analysis.sig)), invert = TRUE)
        # take each path, split into two and unlist them
        pathway.id   <- unlist(lapply(1:length(m),function(n) split(m[[n]][1], " ")))
        pathway.name <- unlist(lapply(1:length(m),function(n) split(m[[n]][2], " ")))
        rownames(analysis.stats) <- pathway.id
        
        analysis.results<- analysis$greater
        
        # Remove gene sets without zero enrichments
        analysis.results <- analysis.results[complete.cases(analysis.results), ]
        
        # Extract Pathway ID and Names
        m <- regmatches(rownames(analysis.results), regexpr(" ", rownames(analysis.results)), invert = TRUE)
        # take each path, split into two and unlist them
        pathway.id   <- unlist(lapply(1:length(m), function(n) split(m[[n]][1], " ")))
        pathway.name <- unlist(lapply(1:length(m), function(n) split(m[[n]][2], " ")))
        
        # Create top table
        colnames(analysis.results)
        toptable <- data.frame(pathway.id, pathway.name, analysis.results[,1:5])
        toptable <- toptable[order(toptable$p.val),]
        rownames(toptable) <- NULL
        colnames(toptable) <- NULL
        
        # save "Toptable"
        filename <- paste(rundir, "gage_data.json", sep="")
        write(toJSON(list(tops = toptable), digits=I(4)), filename)
        
        # save toptable to a tab delimited file
        colnames(toptable) <- c("PathwayID", "Pathway","PGeomean","StatMean", "PValue","QValue","SetSize") 
        filename <- paste(rundir, "gage_toptable.tsv", sep = "")
        write.table(toptable, filename, col.names=NA, sep = "\t" )
        
        # Creating a heatmap
        get.heatmap(analysis.stats, analysis.type)
        
        filename <- paste(rundir, "gage.RData", sep="")
        save( analysis.type, geo.dataset, analysis,
              Group1, Group1names, Group2,Group2names,
              keggcode.organism,file = filename)
    }else{
        # Exit with error code 1
        print("No Significant Results Found!")
        quit(save = "no", status = 1, runLast = FALSE)
    }
}

#############################################################################
#                        Function Calling                                   #
#############################################################################

if (isdebug) { print("GAGE: GAGE analysis starting...") }
compare.option <- ifelse(comparison.type =="ExpVsCtrl", "unpaired", "paired")

if (comparison.type =="ExpVsCtrl"){
    G2 <- Group2
    G1 <- Group1
} else {
    G2 <- NULL
    G1 <- NULL
}

gage.analysis(dbdata, comparison.type, G2, G1, compare.option)

if (isdebug) { print("GAGE analysis completed!") }
quit(save = "no", status = 0, runLast = FALSE)
