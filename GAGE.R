##Analysis specific to the dengue dataset.

#----------------------Parameters to bare in mind-----------------------------

# Function to include the following arguments: 
# Species= To be fed into bods to get org argument value
# The column that contains the different groups
# The identity of the two groups (have an option to merge groups)
# The GO dataset to be used.
# The type of gene sets used (kegg.gs is only for humans)



#----------------------Loading the data-------------------------------------


#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite(c("gage","gageData","GO.db", "pathview" ))
library(GEOquery)
library(gage) #Does the analysis
library(gageData) #Lets data be used by gage
library(pathview) #Visualises interaction networks & used to get ENTREZ IDs
library(GO.db) #Loads GO database

#Importing data from GEO
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

#Get dataset with expression info
X <- Table(gse)

#Converting GSE to an expression set object
eset <- GDS2eSet(gse, do.log2=TRUE)


pDat <- pData(eset)





#---------------------------Using the GAGE package------------------------------


data <- getGEO("GDS5093", GSEMatrix = TRUE)
X <- Table(data)
eset <- GDS2eSet(data, do.log2=TRUE)
pDat <- pData(eset)
##Remove probe ID column & convert into data matrix
X1<- X[,-1]
X1_matrix<-data.matrix(X1)
id.map.refseq <- id2eg(ids = X$IDENTIFIER, category = "SYMBOL", org = "hsa")
for (i in 1:length(id.map.refseq[,1])){
  if (id.map.refseq[i,1] == X1_matrix[i,1]){
    X1_matrix[i,1]<-id.map.refseq[i,2]
  }
}
X1_matrix<-X1_matrix[complete.cases(X1_matrix),]
##Make first column rownames
GEOdataset <- X1_matrix[,-1]
rownames(GEOdataset) <- X1_matrix[,1]
##Convert to numerical matrix (for gage function)
class(GEOdataset) <- "numeric"  
#Get positions of specific, in vector form
cn=colnames(GEOdataset)
##Find out which group each sample is part of
##Use this to form two groups of samples for further analysis
Group1<-c()
Group2<-c()
for (a in 1:length(pDat$sample)){
  if (pDat$infection[a] == "Dengue virus"){
    Group1<-c(Group1, (grep(pDat$sample[a], cn)))
  }
  if (pDat$infection[a] == "control"){
    Group2<- c(Group2, (grep(pDat$sample[a], cn)))
  }
}

data(kegg.gs)
kg.hsa=kegg.gsets("hsa") #this picks out the human sets
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step
save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object
##Using the gage function to carry out analysis
#name           <- gage(data    ,  genesets used, control group, experimental group)
GEOdatasetA.kegg.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, same.dir = F, compare = 'unpaired')
GEOdatasetB.kegg.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group1, samp = Group2, same.dir = F, compare = 'unpaired')

GEOdatasetA.kegg.p.sig<-sigGeneSet(GEOdatasetA.kegg.p, outname="GEOdatasetA.kegg")
GEOdatasetB.kegg.p.sig<-sigGeneSet(GEOdatasetB.kegg.p, outname="GEOdatasetB.kegg")

GEOdatasetA.kegg.2d.sig<-GEOdatasetA.kegg.p.sig[,grep("^stats.GSM", names(GEOdatasetA.kegg.p.sig), value=TRUE)]

heatmap.2(as.matrix(GEOdatasetA.kegg.2d.sig[1:20,]), dendrogram = "none", key=T, keysize=1.5, main = "Top 20 Enriched Gene Sets", trace="none", density.info="none")

GEOdatasetB.kegg.2d.sig<-GEOdatasetB.kegg.p.sig[,grep("^stats.GSM", names(GEOdatasetB.kegg.p.sig), value=TRUE)]

heatmap.2(as.matrix(GEOdatasetB.kegg.2d.sig[1:20,]), dendrogram = "none", key=T, keysize=1.5, main = "Top 20 Enriched Gene Sets", trace="none", density.info="none")





#------------------------Data Preparation----------------------------------------


##Creating table of organisms IDs
data(bods)
bods<-as.data.frame(bods, stringsAsFactors= TRUE )
latin_names<-c("Anopheles","Arabidopsis thaliana", "Bos taurus", "Caenorhabditis elegans", "Canis lupus familiaris", "Drosophila melanogaster", "Danio rerio", "E coli", "Escherichia coli O157", "Gallus gallus", "Homo sapiens", "Mus musculus", "Macaca mulatta", "Anopheles gambiae", "Pan", "Rattus norvegicus", "Saccharomyces cerevisiae", "Sus scrofa", "Xenopus laevis	") 
bods2<-cbind(bods, latin_names)



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
cn=colnames(GEOdataset)
Group1<-c()
Group2<-c()

for (a in 1:length(pDat$sample)){
  if (pDat$infection[a] == "Dengue virus"){
    Group1<-c(Group1, (grep(pDat$sample[a], cn)))
  }
  if (pDat$infection[a] == "control"){
    Group2<- c(Group2, (grep(pDat$sample[a], cn)))
  }
}

data(kegg.gs)
kg.hsa=kegg.gsets("hsa") #this picks out the human sets
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step
save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object
##Using the gage function to carry out analysis
#name           <- gage(data    ,  genesets used, control group, experimental group)
GEOdataset.kegg.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, compare= 'unpaired')






 
#-------------------Generally Applicable Gene-set Enrichment (GAGE)--------------------


#Get positions of specific, in vector form
cn=colnames(GEOdataset)

##Find out which group each sample is part of
##Use this to form two groups of samples for further analysis


Group1<-c()
Group2<-c()

for (a in 1:length(pDat$sample)){
  if (pDat$infection[a] == "Dengue virus"){
    Group1<-c(Group1, (grep(pDat$sample[a], cn)))
  }
  if (pDat$infection[a] == "control"){
    Group2<- c(Group2, (grep(pDat$sample[a], cn)))
  }
}



##Load gene sets
data(kegg.gs)


##Loading human kegg gene sets

kg.hsa=kegg.gsets("hsa") #this picks out the human sets

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step

save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object


##Using the gage function to carry out two-way analysis

GEOdataset.kegg.2d.p <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, same.dir = F, compare='unpaired')




#---------------------Visualisation & Results-------------------------------------

#Table for two-analysis (all gene sets)
write.table(GEOdataset.kegg.2d.p$greater, file = "Test1", sep = "\t")



#Table show top significant gene sets (for 2 way analysis)
write.table(GEOdataset.kegg.2d.sig$greater, file = "Test2", sep = "\t")



##Producing line summaries

#Returns number of up and down regulated gene sets
GEOdataset.kegg.sig<-sigGeneSet(GEOdataset.kegg.p, outname="GEOdataset.kegg")


##Returns number of two-direction enriched gene sets
GEOdataset.kegg.2d.sig<-sigGeneSet(GEOdataset.kegg.2d.p, outname="GEOdataset.kegg")


##Producing the heatmap
GEOdataset.kegg.2d.sig<-GEOdataset.kegg.2d.sig[,grep("^stats.GSM", names(GEOdataset.kegg.2d.sig), value=TRUE)]
#install.packages("gplots")
library(gplots)
heatmap.2(as.matrix(GEOdataset.kegg.2d.sig[1:20,]), dendrogram = "none", key=T, keysize=1.5, main = "Top 20 Enriched Gene Sets", trace="none", density.info="none")



##Interaction networks

GEOdataset.d<-GEOdataset[, Group1]-rowMeans(GEOdataset[,Group2])

##For upregulated gene pathways
sel <- GEOdataset.kegg.p$greater[, "q.val"] < 0.1 & !is.na(GEOdataset.kegg.p$greater[, "q.val"])
path.ids <- rownames(GEOdataset.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8) 

##Produces  top 3 interaction networks
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))


##For down regulated gene pathways
sel2 <- GEOdataset.kegg.p$less[, "q.val"] < 0.1 & !is.na(GEOdataset.kegg.p$greater[, "q.val"])
path.ids3 <- rownames(GEOdataset.kegg.p$greater)[sel]
path.ids4 <- substr(path.ids, 1, 8) 

##Produces  top 3 interaction networks
pv.out.list2 <- sapply(path.ids4[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))











