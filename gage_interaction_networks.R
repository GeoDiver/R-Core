#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript gage_interaction_networks.R --dbrdata ~/Desktop/kegg.RData --rundir ~/Desktop/ --pathid "hsa00480"
# ---------------------------------------------------------#

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
suppressMessages(library("pathview"))      # Interaction networks & used to get ENTREZ IDs

#-------------------------------Set parsers---------------------------------------

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Paeameters
parser <- add_argument(parser, "--dbrdata",
                       help = "Downloaded GEO dataset full path")
parser <- add_argument(parser, "--rundir",
                       help = "The output directory where graphs get saved")
parser <- add_argument(parser, "--pathid",
                       help = "Interaction Network Path ID")

# allows arguments to be run via the command line
argv <- parse_args(parser)

#############################################################################
#                        Command Line Arguments Retrieval                   #
#############################################################################

# General Parameters
rundir          <- argv$rundir
dbrdata         <- argv$dbrdata
pathid          <- argv$pathid

#############################################################################
#                          Loading Saved Dataset                            #
#############################################################################

filename <- argv$dbrdata #paste(rundir,"kegg.RData", sep = "")

if (file.exists(filename)){
    load(file = filename)
}else{
    print("ERROR:File not found")
    q(save = "no")
}    


#############################################################################
#                          Interaction Networks                             #
#############################################################################

if(analysis.type =="ExpVsCtrl"){
    
    #Find expression change between experimental group and control
    GEOdataset.diff<-GEOdataset[, Group1] - rowMeans(GEOdataset[,Group2])
    
#    sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
#     path.ids <- rownames(analysis$greater)[sel]
#     path.ids2 <- substr(path.ids, 1, 8) 

    pathview(gene.data = GEOdataset.diff[,1:2], pathway.id = pathid, 
             species = keggcode.organism, out.suffix = "gage_pathway")
    
}
if(analysis.type =="ExpVsExp"){
    
#     sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
#     path.ids <- rownames(analysis$greater)[sel]
#     path.ids2 <- substr(path.ids, 1, 8) 
    
    ##Interaction pathways for experimental group 1
    pathview(gene.data = GEOdataset[,Group1names][,1:2], pathway.id = pathid, 
             species = keggcode.organism, out.suffix = "gage_pathway")
    
    ##Interaction pathways for experimental group 2
    pathview(gene.data = GEOdataset[,Group2names][,1:2], pathway.id = pathid, 
             species = keggcode.organism, out.suffix = "gage_pathway")
    
}


