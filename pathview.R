#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : pathview.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : create interaction network    #
# Rscript Pathview.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.rData --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --rundir "~/Desktop/" --dev TRUE
# ---------------------------------------------------------#

# Load dependencies

library("argparser")    # Argument passing
library("pathview")     # Creates Heatmap



# Add parsers

parser <- add_argument(parser, "--gagedata", 
                       help="The GAGE results full path")

parser <- add_argument(parser, "--analysis.type",
                       help="Genesets used and whether or not B should be treated as a control")

parser <- add_argument(parser, "--geneset", 
                       help="The KEGG geneset code")

parser <- add_argument(parser, "--rundir", 
                       help="The output directory where graphs get saved")




parser <- add_argument(parser, "--popA",
                       help = "GroupA - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--popB",
                       help = "GroupB - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--factor", 
                       help="Factor type to be classified by")



# allows arguments to be run via the command line
argv <- parse_args(parser)


#General parameters
gagedata   <- argv$gagedata
analysis.type <- argv$analysis.type
factor.type <- argv$factor 
population1     <- unlist(strsplit(argv$popA, ","))
population2     <- unlist(strsplit(argv$popB, ","))



#Loading and preparing data
data <- load(gagedata)







if(analysis.type =="KEGG_ExpVsCtrl"){
  #Find expression change between experimental group and control
  GEOdataset.d<-GEOdataset[, Group1] - rowMeans(GEOdataset[,Group2])
  
  ########## COMMON ##########
  sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
  path.ids <- rownames(analysis$greater)[sel]
  path.ids2 <- substr(path.ids, 1, 8) 
  ########## COMMON ##########
  
  ##Produces  top 3 interaction networks (from 2 way analysis)
  pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))
  
}
if(analysis.type =="KEGG_ExpVsExp"){
  
  ########## COMMON ##########
  sel <- analysis$greater[, "q.val"] < 0.1 & !is.na(analysis$greater[, "q.val"])
  path.ids <- rownames(analysis$greater)[sel]
  path.ids2 <- substr(path.ids, 1, 8) 
  ########## COMMON ##########
  
  ##Interaction pathways for experimental group 1
  pv.out.list2 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset[,Group1names][,1:2], pathway.id = pid, species = "hsa"))
  
  ##Interaction pathways for experimental group 2
  pv.out.list3 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset[,Group2names][,1:2], pathway.id = pid, species = "hsa"))
  
}



