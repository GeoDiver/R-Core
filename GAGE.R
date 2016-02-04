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
#biocLite(c("GEOquery","gage","gageData","GO.db", "pathview","pheatmap","RColorBrewer"))
library(GEOquery)
library(gage)         #Does the analysis
library(gageData)     #Lets data be used by gage
library(pathview)     #Visualises interaction networks & used to get ENTREZ IDs
library(GO.db)        #Loads GO database
library(pheatmap)     #Used to create heatmap
library(RColorBrewer) #Color palette for heatmap


#Importing data from GEO
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

#Get dataset with expression info
X <- Table(gse)

#Converting GSE to an expression set object
eset <- GDS2eSet(gse, do.log2=TRUE)


pDat <- pData(eset)


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




#---------------------GAGE analysis for group 1-------------------------------------

##Using the gage function to carry out two-way analysis

keggresults_group1 <- gage(GEOdataset, gsets = kegg.gs, ref = Group2, samp = Group1, same.dir = F, compare='unpaired')


##Producing line summaries


##Returns number of two-direction significantly enriched gene sets
keggresults_group1_sig<-sigGeneSet(keggresults_group1)

##Formatting and preparation for heatmap

keggresults_group1_sig<-as.data.frame(keggresults_group1_sig)
keggresults_group1_stats<-keggresults_group1_sig[,grep("^stats.GSM", names(keggresults_group1_sig), value=TRUE)]



##Interaction networks

GEOdataset.d<-GEOdataset[, Group1]-rowMeans(GEOdataset[,Group2])

##For upregulated gene pathways
sel <- keggresults_group1$greater[, "q.val"] < 0.1 & !is.na(keggresults_group1$greater[, "q.val"])
path.ids <- rownames(keggresults_group1$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8) 

##Produces  top 3 interaction networks (from 2 way analysis)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d[,1:2], pathway.id = pid, species = "hsa"))

##Results table

Group1_results<-keggresults_group1$greater

##Remove gene sets without zero enrichments
Group1_results<-Group1_results[complete.cases(Group1_results),]


#---------------------GAGE analysis for group 2-------------------------------------

##Using the gage function to carry out two-way analysis

keggresults_group2 <- gage(GEOdataset, gsets = kegg.gs, ref = Group1, samp = Group2, same.dir = F, compare='unpaired')

##Returns number of two-direction significantly enriched gene sets
keggresults_group2_sig<-sigGeneSet(keggresults_group2)


##Formatting and preparation for heatmap

keggresults_group2_sig<-as.data.frame(keggresults_group2_sig)
keggresults_group2_stats<-keggresults_group2_sig[,grep("^stats.GSM", names(keggresults_group2_sig), value=TRUE)]


##Interaction networks

GEOdataset.d2<-GEOdataset[, Group2]-rowMeans(GEOdataset[,Group1])

##For upregulated gene pathways
sel2 <- keggresults_group2$greater[, "q.val"] < 0.1 & !is.na(keggresults_group1$greater[, "q.val"])
path.ids3 <- rownames(keggresults_group1$greater)[sel]
path.ids4 <- substr(path.ids, 1, 8) 

##Produces  top 3 interaction networks (from 2 way analysis)
pv.out.list2 <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = GEOdataset.d2[,1:2], pathway.id = pid, species = "hsa"))


##Results table

Group2_results<-keggresults_group2$greater

##Remove gene sets without zero enrichments
Group2_results<-Group2_results[complete.cases(Group2_results),]




#---------------------Making the heatmap-------------------------------------

##Combining tables 

allsamples<-merge(keggresults_group1_stats,keggresults_group2_stats, by= "row.names", all=FALSE)
allsamples2 <- allsamples[,-1]
rownames(allsamples2) <- allsamples[,1]


##Creating a heatmap

allsamples2<-t(allsamples2)
row.names(allsamples2)<-gsub("(stats.)", "", row.names(allsamples2))
col.pal <- RColorBrewer::brewer.pal(9, "Reds")
annotation_col <- data.frame( Infection = pDat[,2])
rownames(annotation_col) = pDat[,1]

pheatmap::pheatmap(t(allsamples2), 
                   cluster_row = T,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   color = col.pal, 
                   fontsize = 6.5,
                   fontsize_row=6, 
                   fontsize_col = 6,
                   gaps_col=length(Group1))





##make a named list or dataframe (for mean values)