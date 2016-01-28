#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : GAGE.R 								   #
# Authors       : Marian, Anissa 						   #
# Description   : Generally Applicable Gene Set 		   # 
# ---------------------------------------------------------#


#### TO RUN IN CLI (modify for GAGE) ####
# Rscript GAGE.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005
####


#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")
# biocLite(c("gage","gageData","GO.db", "pathview" ))
# install.packages("argparser")
# install.packages('Cairo')
# install.packages('GEOquery')

### note: add any packages and sort the libraries alphabetically ###

library('argparser')		# Allows the Rscript to run with CLI arguments
library('Cairo') 			# Allows storing plots as .png files
library('gage') 			# Does the analysis
library('gageData') 		# Lets data be used by gage
library('GEOquery')			# 
library('GO.db') 			# Downloads GO datasets
library('pathview') 		# Visualises interaction networks & used to get ENTREZ IDs


#############################################################################
#                        Command Line Arguments                             #
#############################################################################

## Currently: Analysis specific to the dengue dataset ##

#----------------------Parameters to bear in mind-----------------------------#

# Function to include the following arguments: 
# Species= To be fed into bods to get org argument value
# The column that contains the different groups
# The identity of the two groups (have an option to merge groups)
# The GO dataset to be used.
# The type of gene sets used (kegg.gs is only for humans)
# paired or unpaired


#-----------------------------Parsers----------------------------------------#

# # set parsers for all input arguments
# parser <- arg_parser("This parser contains the input arguments")
# parser <- add_argument(parser, "--accession"    , help="input file")    # GEO Accession ID
# parser <- add_argument(parser, "--factor"     , help="input file")    # Factor type to be classified by
# parser <- add_argument(parser, "--popA", nargs='+', help="input file")    # GroupA - all the selected phenotypes (atleast one)
# parser <- add_argument(parser, "--popB", nargs='+', help="input file")    # GroupB - all the selected phenotypes (atleast one)
# parser <- add_argument(parser, "--popname1"     , help="input file")    # name for GroupA
# parser <- add_argument(parser, "--popname2"     , help="input file")    # name for GroupB
# parser <- add_argument(parser, "--topgenecount"   , help="input file")    # number of top genes to be used
# parser <- add_argument(parser, "--foldchange"   , help="input file")    # fold change cut off
# parser <- add_argument(parser, "--thresholdvalue" , help="input file")    # threshold value cut off
# parser <- add_argument(parser, "--working_dir"    , help="input file")    # GEO Accession ID

# # allow arguments to be run via the command line
# argv <- parse_args(parser)

# # output directory (CHANGE)
# output.dir     <- argv$working_dir



#----------------------Loading the data-------------------------------------#

# Importing data from GEO
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

# Get dataset with expression info
X <- Table(gse)

# Converting GSE to an expression set object
eset <- GDS2eSet(gse, do.log2=TRUE)

pDat <- pData(eset)


#############################################################################
#                      Using the GAGE package                   	        #
#############################################################################


#------------------------Data Preparation-----------------------------------#

##Remove probe ID column & convert into data matrix
X1<- X
X1<- X[,-1]
X1_matrix<-data.matrix(X1)


##Create two column table containing entrez IDs for geodataset
id.map.refseq <- id2eg(ids = X$IDENTIFIER, category = "SYMBOL", org = "hsa")
#data(bods) - contains values  for 'org' argument. 


##Replace gene symbols with ENTREZ ID in dataset matrix
for (i in 1:length(id.map.refseq[,1])){
  if (id.map.refseq[i,1] == X1_matrix[i,1]){
    X1_matrix[i,1]<-id.map.refseq[i,2]
  }
}

##Remove rows without ENTREZ IDs
X1_matrix<-X1_matrix[complete.cases(X1_matrix),]

##Make first column rownames
GEOdataset <- X1_matrix[,-1]
rownames(GEOdataset) <- X1_matrix[,1]
##Convert to numerical matrix (for gage function)
class(GEOdataset) <- "numeric"  

 
#-------------------Generally Applicable Gene-set Enrichment (GAGE)--------------------


#Get positions of specific, in vector form
cn=colnames(X1_matrix2)

##Find out which group each sample is part of
##Use this to form two groups of samples for further analysis


Group1<-c()
Group2<-c()

for (a in 1:length(pDat$sample)){
  if (pDat$infection[a] == "Dengue virus"){
    Group1_pos<- grep(pDat$sample[a], cn)
    Group1<-c(Group1, Group1_pos)
  }
  if (pDat$infection[a] == "control"){
    Group2<- c(Group2, (grep(pDat$sample[a], cn)))
  }
}

##Load gene sets
data(kegg.gs)
data(go.gs)


##Loading human kegg gene sets

kg.hsa=kegg.gsets("hsa") #this picks out the human sets

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step

save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object

##Loading human GO gene sets

go.hs=go.gsets(species="human")
go.bp=go.hs$go.sets[go.hs$go.subs$BP]
go.mf=go.hs$go.sets[go.hs$go.subs$MF]
go.cc=go.hs$go.sets[go.hs$go.subs$CC]
save(go.bp, go.mf, go.cc, file="go.hs.gsets.RData")
#for Bioconductor species supported by go.gsets function:
#data(bods)
#print(bods)



##Using the gage function to carry out analysis


#name           <- gage(data    ,  genesets used, control group, experimental group)
GEOdataset.kegg.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, compare= 'unpaired')

GEOdataset.go.p <- gage(GEOdataset, gsets = go.gs, ref = Group2, samp = Group1,compare= 'unpaired')

GEOdataset.bp.p <- gage(GEOdataset, gsets = go.bp,ref = Group2, samp = Group1, compare= 'unpaired')




#---------------------Visualisation & Results-------------------------------------

##Producing tables

#Carry out two-analysis (only for Kegg gene sets)
GEOdataset.kegg.2d.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, same.dir = F, compare='unpaired')

#Table for two-analysis
write.table(GEOdataset.kegg.2d.p$greater, file = "GEOdataset.kegg.2d.p.txt", sep = "\t")

#Table for one-way analysis (upregulated and downgregulated gene sets)
write.table(rbind(GEOdataset.kegg.p$greater, GEOdataset.kegg.p$less), file = "GEOdataset.kegg.p.txt", sep = "\t")

#Table show top significant gene sets (for 2 way analysis)
write.table(GEOdataset.kegg.2d.sig$greater, file = "GEOdataset.kegg.2d.sig.txt", sep = "\t")

#For 1 way analysis
write.table(rbind(GEOdataset.kegg.sig$greater, GEOdataset.kegg.sig$less), file = "GEOdataset.kegg.sig.txt", sep = "\t")



##Producing line summaries

#Returns number of up and down regulated gene sets
GEOdataset.kegg.sig<-sigGeneSet(GEOdataset.kegg.p, outname="GEOdataset.kegg")

##Its heatmap
sigGeneSet(GEOdataset.kegg.p, outname="GEOdatasetUP.kegg", heatmap= TRUE)

##Returns number of two-direction enriched gene sets
GEOdataset.kegg.2d.sig<-sigGeneSet(GEOdataset.kegg.2d.p, outname="GEOdataset.kegg")




## Creating a heatmap  #Need to convert back to gene symbol
gs=unique(unlist(kegg.gs[rownames(GEOdataset.kegg.p$greater)[1:3]]))


#Extract data for essential member genes in a gene set
#Creates a table, and the row headers consist of the samples
#In the order Group2, Group1
essData=essGene(gs, GEOdataset, ref =Group2, samp =Group1, compare="unpaired")


#Provide the new positions of the sample groups
ref1=1:length(Group2)

samp1=(length(ref1)+1): length(essData[1,])

#Produces heatmaps for top 3 up-regulated KEGG pathways in a batch

for (gs in rownames(GEOdataset.kegg.p$greater)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
           samp = samp1, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}



#Produces heatmaps for top 3 down-regulated KEGG pathways in a batch
#Too many genes!!!!! Consider using a cut off?

for (gs in rownames(GEOdataset.kegg.p$less)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = GEOdataset, ref = Group2,
           samp = Group1, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


##Non-redundant gene sets (can feasibly do table and heatmap??)

GEOdataset.kegg.esg.up <- esset.grp(GEOdataset.kegg.p$greater, GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1,test4up = T, output = T, outname = "GEOdataset.kegg.up", make.plot = F, compare="unpaired")

GEOdataset.kegg.esg.dn <- esset.grp(GEOdataset.kegg.p$less, GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, test4up = F, output = T, outname = "GEOdataset.kegg.dn", make.plot = F, compare="unpaired")


##Interaction networks

GEOdataset.d<-GEOdataset[, Group1]-rowMeans(GEOdataset[,Group2])

##For upregulated gene pathways
sel <- GEOdataset.kegg.p$greater[, "q.val"] < 0.1 & !is.na(GEOdataset.kegg.p$greater[, "q.val"])
path.ids <- rownames(GEOdataset.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8) 
##Produces  top 10 interaction networks
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))

##For down regulated gene pathways
sel2 <- GEOdataset.kegg.p$less[, "q.val"] < 0.1 & !is.na(GEOdataset.kegg.p$greater[, "q.val"])
path.ids3 <- rownames(GEOdataset.kegg.p$greater)[sel]
path.ids4 <- substr(path.ids, 1, 8) 

##Produces  top 10 interaction networks
pv.out.list2 <- sapply(path.ids4[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))

