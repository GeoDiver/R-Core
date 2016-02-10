#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
    # Rscript GAGE.R --accession GDS5093 --factor "infection" --outputdir "~/Desktop/" --organism "hsa"
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

# silence library loading messages on command line
suppressMessages(library("GEOquery"))
suppressMessages(library("gage"))
suppressMessages(library("gageData"))
suppressMessages(library("pathview"))
suppressMessages(library("GO.db"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("pheatmap"))

# load required libraries
library(argparser)    # Argument passing
library(gage)         # Does the analysis
library(gageData)     # Lets data be used by gage
library(GEOquery)     # GEO dataset Retrieval
library(GO.db)        # Loads GO database
library(pathview)     # Visualises interaction networks & used to get ENTREZ IDs
library(pheatmap)     # Used to create heatmap
library(RColorBrewer) # Color palette for heatmap

#-------------------------------Set parsers---------------------------------------

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

parser <- add_argument(parser, "--accession", 
                       help="Accession Number of the GEO Database")
parser <- add_argument(parser, "--factor", 
                       help="Factor type to be classified by")
parser <- add_argument(parser, "--organism", 
                       help="Organism which the data comes from")
parser <- add_argument(parser, "--outputdir", 
                       help="The output directory where graphs get saved")

# allows arguments to be run via the command line
argv <- parse_args(parser)


#------------------------------Set Parameters-------------------------------------

# General Parameters
output.dir  <- argv$outputdir

# Sample Parameters
accession   <- argv$accession # GDS5093
factor   <- argv$factor    # "factor"
pop.colour1 <- "#b71c1c"      # Red
pop.colour2 <- "#0d47a1"      # Blue
organism    <- argv$organism  # "hsa"

rundir      <- "/Users/sureshhewapathirana/Desktop/"
accession   <- "GDS5093"
factor      <- "infection"
organism    <- "hsa"
#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################

# Importing data from GEO
gse  <- getGEO(accession, GSEMatrix = TRUE)

# Converting GSE to an expression set object
eset <- GDS2eSet(gse, do.log2=TRUE)

# Get dataset with expression info
Y    <- Table(gse)

pDat <- pData(eset)

# Annotation column for heatmap and grouping
annotation_col <- pDat[factor]
rownames(annotation_col) = pDat[,1]



#------------------------Data Preparation----------------------------------------

## Creating table of organisms IDs
data(bods)

bods        <- as.data.frame(bods, stringsAsFactors= TRUE )
latin_names <- c("Anopheles","Arabidopsis thaliana", "Bos taurus", "Caenorhabditis elegans", "Canis lupus familiaris", "Drosophila melanogaster", "Danio rerio", "E coli", "Escherichia coli O157", "Gallus gallus", "Homo sapiens", "Mus musculus", "Macaca mulatta", "Anopheles gambiae", "Pan", "Rattus norvegicus", "Saccharomyces cerevisiae", "Sus scrofa", "Xenopus laevis	") 
bods2       <- cbind(bods, latin_names)


## Remove probe ID column & convert into data matrix
Y1        <- Y
Y1        <- Y[,-1]
Y1_matrix <-data.matrix(Y1)


## Create two column table containing entrez IDs for geodataset
id.map.refseq <- id2eg(ids = Y$IDENTIFIER, category = "SYMBOL", org = organism)
#data(bods) - contains values  for 'org' argument. 


##Replace gene symbols with ENTREZ ID in dataset matrix
X1_matrix[,1]<- id.map.refseq[,2]

## Remove rows without ENTREZ IDs
Y1_matrix<-Y1_matrix[complete.cases(Y1_matrix),]

## Make first column rownames
GEOdataset <- Y1_matrix[,-1]
rownames(GEOdataset) <- Y1_matrix[,1]

## Convert to numerical matrix (for gage function)
class(GEOdataset) <- "numeric"  


## Get group position and sample names
Group1<- which(annotation_col$infection == "Dengue virus")
Group1names<-rownames(annotation_col)[annotation_col$infection == "Dengue virus" ]
Group2<- which(annotation_col$infection == "control")
Group2names<-rownames(annotation_col)[annotation_col$infection == "control" ]

##Loading kegg sets

data(kegg.gs)
kg.hsa=kegg.gsets(organism) #this picks out the human sets
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step
save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object



##Loading GO sets

#BP = Biological Process MF = molecular function CC = cellular component
go.hs=go.gsets(species="human")  #use species column of bods2
go.bp=go.hs$go.sets[go.hs$go.subs$BP]
go.mf=go.hs$go.sets[go.hs$go.subs$MF]
go.cc=go.hs$go.sets[go.hs$go.subs$CC]
save(go.bp, go.mf, go.cc, file="go.hs.gsets.RData")



#############################################################################
#               GAGE analysis for experimental vs  control                  #
#############################################################################


##Kegg gene sets
#----------------


##Using the gage function to carry out two-way analysis

keggresults_analysis1 <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, same.dir = F, compare='unpaired')


##Returns number of two-direction significantly enriched gene sets

keggresults_analysis1_sig<-sigGeneSet(keggresults_analysis1)


##Formatting and preparation for heatmap

keggresults_analysis1_sig<-as.data.frame(keggresults_analysis1_sig)
keggresults_analysis1_stats<-keggresults_analysis1_sig[,grep("^stats.GSM", names(keggresults_analysis1_sig), value=TRUE)]



##Interaction networks


#Find expression change between experimental group and control
GEOdataset.d<-GEOdataset[, Group1]-rowMeans(GEOdataset[,Group2])


sel <- keggresults_analysis1$greater[, "q.val"] < 0.1 & !is.na(keggresults_analysis1$greater[, "q.val"])
path.ids <- rownames(keggresults_analysis1$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8) 

##Produces  top 3 interaction networks (from 2 way analysis)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))


##Results table

Analysis1_results<-keggresults_analysis1$greater

##Remove gene sets without zero enrichments
Analysis1_results<-Analysis1_results[complete.cases(Analysis1_results),]



##Creating a heatmap

Analysis1_heatmap<-t(keggresults_analysis1_stats)
Analysis1_heatmap<-Analysis1_heatmap[,1:20]
row.names(Analysis1_heatmap)<-gsub("(stats.)", "", row.names(Analysis1_heatmap))
col.pal <- RColorBrewer::brewer.pal(9, "Reds")


pheatmap::pheatmap(t(Analysis1_heatmap), 
                   cluster_row = F,
                   cluster_cols = T,
                   color = col.pal, 
                   fontsize = 6.5,
                   fontsize_row=6, 
                   fontsize_col = 6)



#Gene ontology sets
#------------------


#arguments: go.cc, go.mf, go.bp

GO_ExpVsCtrl <- function(set_type){
  
  
  ##Using the gage function to carry out two-way analysis
  
  GOresults_analysis1 <- gage(GEOdataset, gsets = set_type, ref = Group2, samp = Group1, same.dir = F, compare='unpaired')
  
  
  ##Returns number of two-direction significantly enriched gene sets
  GOresults_analysis1_sig<-sigGeneSet(GOresults_analysis1)
  
  
  ##Formatting and preparation for heatmap
  GOresults_analysis1_sig<-as.data.frame(GOresults_analysis1_sig)
  GOresults_analysis1_stats<-GOresults_analysis1_sig[,grep("^stats.GSM", names(GOresults_analysis1_sig), value=TRUE)]
  
  ##Results table
  Analysis1_results<-GOresults_analysis1$greater
  
  ##Remove gene sets without zero enrichments
  Analysis1_results<-Analysis1_results[complete.cases(Analysis1_results),]
  
  
  
  ##Creating a heatmap
  
  Analysis1_heatmap<-t(GOresults_analysis1_stats)
  Analysis1_heatmap<-Analysis1_heatmap[,1:20]
  row.names(Analysis1_heatmap)<-gsub("(stats.)", "", row.names(Analysis1_heatmap))
  col.pal <- RColorBrewer::brewer.pal(9, "Reds")
  
  
  pheatmap::pheatmap(t(Analysis1_heatmap), 
                     cluster_row = F,
                     cluster_cols = T,
                     color = col.pal, 
                     fontsize = 6.5,
                     fontsize_row=6, 
                     fontsize_col = 6)
}




#############################################################################
#          GAGE analysis for two experimental groups                        #
#############################################################################


##Kegg gene sets
#---------------


##Using the gage function to carry out two-way analysis

keggresults_analysis2 <- gage(GEOdataset, gsets= kegg.gs, ref=NULL , samp=NULL, same.dir = F)

##Returns number of two-direction significantly enriched gene sets

keggresults_analysis2_sig<-sigGeneSet(keggresults_analysis2)


##Formatting and preparation for heatmap

keggresults_analysis2_sig<-as.data.frame(keggresults_analysis2_sig)
keggresults_analysis2_stats<-keggresults_analysis2_sig[,grep("^stats.GSM", names(keggresults_analysis2_sig), value=TRUE)]


##Interaction networks

#subset GEOdataset with experimental group 1
Exp1<- GEOdataset[,Group1names]


#subset GEOdataset with experimental group 2
Exp2<- GEOdataset[,Group2names]


##Interaction pathways 
sel2 <- keggresults_analysis2$greater[, "q.val"] < 0.1 & !is.na(keggresults_analysis2$greater[, "q.val"])
path.ids3 <- rownames(keggresults_analysis2$greater)[sel]
path.ids4 <- substr(path.ids, 1, 8) 

##Interaction pathways for experimental group 1
pv.out.list2 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = Exp1[,1:2], pathway.id = pid, species = "hsa"))


##Interaction pathways for experimental group 2

pv.out.list3 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = Exp2[,1:2], pathway.id = pid, species = "hsa"))


##Results table

Analysis2_results<-keggresults_analysis2$greater

##Remove gene sets without zero enrichments
Analysis2_results<-Analysis2_results[complete.cases(Analysis2_results),]

##Creating a heatmap

Analysis2_heatmap<-t(keggresults_analysis2_stats)
Analysis2_heatmap<-Analysis2_heatmap[,1:20]
row.names(Analysis2_heatmap)<-gsub("(stats.)", "", row.names(Analysis2_heatmap))
col.pal <- RColorBrewer::brewer.pal(9, "Reds")

pheatmap::pheatmap(t(Analysis2_heatmap), 
                   cluster_row = F,
                   cluster_cols = T,
                   annotation_col = annotation_col,
                   color = col.pal, 
                   fontsize = 6.5,
                   fontsize_row=6, 
                   fontsize_col = 6,
                   gaps_col=length(Group1))




##Gene ontology sets
#-------------------

#arguments: go.cc, go.mf, go.bp

GO_ExpVsExp <- function(set_type){
  ##Using the gage function to carry out two-way analysis
  
  GOresults_analysis2 <- gage(GEOdataset, gsets= set_type, ref=NULL , samp=NULL, same.dir = F)
  
  ##Returns number of two-direction significantly enriched gene sets
  
  GOresults_analysis2_sig<-sigGeneSet(GOresults_analysis2)
  
  
  ##Formatting and preparation for heatmap
  
  GOresults_analysis2_sig<-as.data.frame(GOresults_analysis2_sig)
  GOresults_analysis2_stats<-GOresults_analysis2_sig[,grep("^stats.GSM", names(GOresults_analysis2_sig), value=TRUE)]
  
  ##Results table
  
  Analysis2_results<-GOresults_analysis2$greater
  
  ##Remove gene sets without zero enrichments
  Analysis2_results<-Analysis2_results[complete.cases(Analysis2_results),]
  
  
  
  ##Creating a heatmap
  
  Analysis2_heatmap<-t(GOresults_analysis2_stats)
  Analysis2_heatmap<-Analysis2_heatmap[,1:20]
  row.names(Analysis2_heatmap)<-gsub("(stats.)", "", row.names(Analysis2_heatmap))
  col.pal <- RColorBrewer::brewer.pal(9, "Reds")
  
  
  
  pheatmap::pheatmap(t(Analysis2_heatmap), 
                     cluster_row = F,
                     cluster_cols = T,
                     annotation_col = annotation_col,
                     color = col.pal, 
                     fontsize = 6.5,
                     fontsize_row=6, 
                     fontsize_col = 6,
                     gaps_col=length(Group1))
}