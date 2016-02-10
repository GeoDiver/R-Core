#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript gage.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.RData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control"  --comparisontype ExpVsCtrl --genesettype KEGG --geotype BP --dev TRUE
# ---------------------------------------------------------#

## Analysis specific to the dengue dataset

#----------------------Parameters to bare in mind-----------------------------

# Function to include the following arguments: 
# Species= To be fed into bods to get org argument value
# The column that contains the different groups
# The identity of the two groups (have an option to merge groups)
# The GO dataset to be used.
# The type of gene sets used (kegg.gs is only for humans)


#############################################################################
#                        Import Libraries                                   #
#############################################################################

# load required libraries and silence library loading messages on command line
suppressMessages(library("argparser"))     # Argument passing
suppressMessages(library("Cairo"))         # Plots saving
suppressMessages(library("gage"))          # Does the analysis
suppressMessages(library("gageData"))      # Lets data be used by gage
suppressMessages(library("GEOquery"))      # GEO dataset Retrieval
suppressMessages(library("GO.db"))         # Loads GO database
suppressMessages(library("jsonlite"))      # Convert R object to JSON format
suppressMessages(library("pathview"))      # Interaction networks & used to get ENTREZ IDs
suppressMessages(library("pheatmap"))      # Used to create heatmap
suppressMessages(library("RColorBrewer"))  # Color palette for heatmap

#-------------------------------Set parsers---------------------------------------

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

# Sample Parameters
parser <- add_argument(parser, "--comparisontype", help = "ExpVsCtrl or ExpVsExp")
parser <- add_argument(parser, "--genesettype",help = "KEGG or GO")
parser <- add_argument(parser, "--geotype",
                       help = "BP - Biological Process or 
                               MF - molecular function 
                               or CC - Cellular Component")

# allows arguments to be run via the command line
argv <- parse_args(parser)

#############################################################################
#                        Command Line Arguments Retrieval                   #
#############################################################################

# General Parameters
rundir          <- argv$rundir
dbrdata         <- argv$dbrdata

if(!is.na(argv$dev)){
	isdebug     <- argv$dev
} else {
	isdebug     <- FALSE
}

# Sample Parameters
accession   <- argv$accession
factor.type <- argv$factor
population1 <- unlist(strsplit(argv$popA, ","))
population2 <- unlist(strsplit(argv$popB, ","))

pop.colour1 <- "#b71c1c"
pop.colour2 <- "#0d47a1"

# Gage parameters
comparison.type <- argv$comparisontype  # "ExpVsCtrl"  # or ExpVsExp
geneset.type    <- argv$genesettype     # "KEGG"  # or "GO"
geo.type        <- argv$geotype         # "BP" # or "MF" or "CC"

#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################
if(isdebug){
	print("GeoDiver is starting")
}

if(isdebug){
	print("Libraries have been loaded")
}

if (file.exists(dbrdata)){
    load(file = dbrdata)
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

if(isdebug){
  print("Datset has been loaded")
  print(paste("Analyzing the factor", factor.type))
  print(paste("for", argv$popA)) 
  print(paste("against", argv$popB))
}

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

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

## Get sample indexes and sample names
Group1<-  which(expression.info[,"population"] == "Group1") 
Group1names<- expression.info[Group1,"Sample"]              
Group2<-  which(expression.info[,"population"] == "Group2") 
Group2names<- expression.info[Group2,"Sample"]  

if(isdebug){
  print("Factors and Populations have been set")
}

#############################################################################
#                            Data Preparation                               #
#############################################################################

## Creating table of organisms IDs
data(bods)
organism    <- as.character(Meta(gse)$sample_organism)

bods        <- as.data.frame(bods, stringsAsFactors= TRUE )
latin_names <- c("Anopheles","Arabidopsis thaliana", "Bos taurus", "Caenorhabditis elegans", 
                 "Canis lupus familiaris", "Drosophila melanogaster", "Danio rerio", "E coli", 
                 "Escherichia coli O157", "Gallus gallus", "Homo sapiens", "Mus musculus", 
                 "Macaca mulatta", "Anopheles gambiae", "Pan", "Rattus norvegicus", 
                 "Saccharomyces cerevisiae", "Sus scrofa", "Xenopus laevis") 
bods2       <- cbind(bods, latin_names)

keggcode.organism <- bods2[which(bods2[,"latin_names"] == organism),"kegg code"]

# Get dataset with expression info
Y           <- Table(gse)

## Remove probe ID column & convert into data matrix
Y1_matrix <-data.matrix(Table(gse)[,-1])

## Create two column table containing entrez IDs for geodataset
id.map.refseq <- id2eg(ids = Y$IDENTIFIER, category = "SYMBOL", org = keggcode.organism)

## Replace gene symbols with ENTREZ ID in dataset matrix
Y1_matrix[,1]<-id.map.refseq[,2]

## Remove rows without ENTREZ IDs
Y1_matrix<-Y1_matrix[complete.cases(Y1_matrix),]

## Make first column rownames
GEOdataset <- Y1_matrix[,-1]
rownames(GEOdataset) <- Y1_matrix[,1]

## Convert to numerical matrix (for gage function)
class(GEOdataset) <- "numeric"  

if(isdebug){
	print("Data Preparation completed")
}

#############################################################################
#                          Gage  Data Loading                               #
#############################################################################

if ('KEGG' == 'KEGG') {
  # Loading kegg sets
  data(kegg.gs)
  kg.hsa  = kegg.gsets(organism)                       #this picks out the human sets
  kegg.gs = kg.hsa$kg.sets[kg.hsa$sigmet.idx]         # no idea but doesn't seem to work without this step
  # filename <- paste(rundir, "kegg.hsa.sigmet.gsets.RData", sep="")
  # save(kegg.gs, file = filename) #saves the human sets as an R object
} else if ('GO' == 'GO') {
  # Loading GO sets
  go.hs = go.gsets(species="human")       # use species column of bods2
  go.bp = go.hs$go.sets[go.hs$go.subs$BP] # BP = Biological Process
  go.mf = go.hs$go.sets[go.hs$go.subs$MF] # MF = molecular function
  go.cc = go.hs$go.sets[go.hs$go.subs$CC] # CC = cellular component
  # filename <- paste(rundir, "go.hs.gsets.RData", sep="")
  # save(go.bp, go.mf, go.cc, file = filename)
}  

if(isdebug){
	print("Gage Data Preparation completed!")
}

#############################################################################
#               Heatmap                  #
#############################################################################

get.heatmap <- function(analysis.stats, heatmap.name){
    
    analysis.heatmap<-t(analysis.stats)
    analysis.heatmap<-analysis.heatmap
    row.names(analysis.heatmap)<-gsub("(stats.)", "", row.names(analysis.heatmap))
    col.pal <- colorRampPalette(rev(
        RColorBrewer::brewer.pal(11, "RdYlGn")))(100)
    
    filename <- paste(rundir, heatmap.name, sep="")

    if(isdebug ){
    	print(paste("Saving heatmap", filename))
    }

    CairoSVG(file = filename)
    pheatmap::pheatmap(t(analysis.heatmap), 
                       cluster_row = F,
                       cluster_cols = T,
                       annotation_col = pclass,
                       color = col.pal, 
                       fontsize = 6.5,
                       fontsize_row=6, 
                       fontsize_col = 6
                       )
    dev.off()

}

#############################################################################
#                        GAGE analysis for KEGG                             #
#############################################################################

kegg.analysis <- function(set.type , analysis.type = "ExpVsCtrl", ref.group = G2, samp.group = G1, compare.option = "paired"){
   
    analysis <- gage(GEOdataset, gsets = kegg.gs,  #set.type, 
                         ref = G2, samp = G1, 
                         same.dir = F, compare = "unpaired" )
    
    # Returns number of two-direction significantly enriched gene sets
    analysis.sig<-sigGeneSet(analysis)
    
    # Formatting and preparation for heatmap
    analysis.sig<-as.data.frame(analysis.sig)
    analysis.stats<-analysis.sig[,grep("^stats.GSM", names(analysis.sig), value=TRUE)]
    
    ##Interaction networks
    
#    if(analysis.type =="ExpVsCtrl"){
#        #Find expression change between experimental group and control
#        GEOdataset.d<-GEOdataset[, Group1] - rowMeans(GEOdataset[,Group2])
#        
#        ########## COMMON ##########
#        sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
#        path.ids <- rownames(analysis$greater)[sel]
#        path.ids2 <- substr(path.ids, 1, 8) 
#        ########## COMMON ##########
#        
#        ##Produces  top 3 interaction networks (from 2 way analysis)
#        pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))
#        
#    }
#    if(analysis.type =="ExpVsExp"){
#        
#        ########## COMMON ##########
#        sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
#        path.ids <- rownames(analysis$greater)[sel]
#        path.ids2 <- substr(path.ids, 1, 8) 
#        ########## COMMON ##########
#        
#        ##Interaction pathways for experimental group 1
#        pv.out.list2 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset[,Group1names][,1:2], pathway.id = pid, species = "hsa"))
#        
#        ##Interaction pathways for experimental group 2
#        pv.out.list3 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset[,Group2names][,1:2], pathway.id = pid, species = "hsa"))
#        
#    }
    
    # Results table
    analysis.results<-analysis$greater
    
    # Remove gene sets without zero enrichments
    analysis.results<-analysis.results[complete.cases(analysis.results),]
    
    # save "Toptable"
    filename <- paste(rundir, "gage_data.json", sep = "")
    write(toJSON(analysis.results, digits=I(4)), filename )
    
    # Creating a heatmap
    get.heatmap(analysis.stats, "gage_heatmap.svg")
}

#############################################################################
#          GAGE analysis for Gene ontology sets                             #
#############################################################################

#arguments: go.cc, go.mf, go.bp

go.analysis <- function(set.type , analysis.type = "ExpVsCtrl", ref.group, samp.group, compare.option = "paired" ){
  
        analysis <- gage(GEOdataset, gsets = set.type, 
                                    ref = ref.group, samp = samp.group, 
                                    same.dir = F, compare= compare.option)
     
    # Returns number of two-direction significantly enriched gene sets
    analysis.sig<-sigGeneSet(analysis)
    
    # Formatting and preparation for heatmap
    analysis.sig <- as.data.frame(analysis.sig)
    analysis.stats<-analysis.sig[,grep("^stats.GSM", names(analysis.sig), value=TRUE)]
    
    # Results table
    analysis.results<-analysis$greater
    
    # Remove gene sets without zero enrichments
    analysis.results<-analysis.results[complete.cases(analysis.results),]
    
    # save "Toptable"
    filename <- paste(rundir, "gage_data.json", sep = "")
    write(toJSON(analysis.results, digits=I(4)), "filename" )
    
    # Creating a heatmap
    get.heatmap(analysis.stats, "gage_heatmap.svg")
}

#############################################################################
#                        Function Calling                                   #
#############################################################################

comp.option <- ifelse(comparison.type =="ExpVsCtrl", "unpaired", "paired")
if(comparison.type =="ExpVsCtrl"){
    G2 <- Group2
    G1 <- Group1
}else{
    G2 <- NULL
    G1 <- NULL
}

if(geneset.type == "KEGG"){
    kegg.analysis(kegg.gs, comparison.type, G2, G1, comp.option)
    if(isdebug ){
    	print("KEGG Analysis completed!")
    }
}

if(geneset.type == "GO"){
    if(geo.type == "BP"){
        go.analysis(go.bp, comparison.type, G2, G1,comp.option)
    } else if(geo.type == "MF"){
        go.analysis(go.mf, comparison.type, G2, G1,comp.option)
    } else if(geo.type == "CC"){
        go.analysis(go.cc, comparison.type, G2, G1,comp.option)
    }
    
    if(isdebug){
    	print("GO Analysis completed!")
	}
}

