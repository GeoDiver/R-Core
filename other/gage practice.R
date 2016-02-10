source("https://bioconductor.org/biocLite.R")


#install.packages('gage')
#biocLite('gage')
#biocLite('gageData')
#biocLite('pathview')
#biocLite('GO.db')

library(GEOquery) #Downloads dataset
library(gage) #Does the analysis
library(pathview) #Visualises interaction networks
library(gageData) #Lets data be used by gage
library(GO.db) #Downloads GO datasets


gse16873 <- getGEO("GSE16873", GSEMatrix = TRUE) #get the data

data(gse16873) #make it into a dataset in R

hn=(1:6)*2-1 #for this example these 2 lines are defining the two groups that will be compared aginst one another
dcis=(1:6)*2

data(kegg.gs) #turns kegg data into a dataset
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis) #groups genes based on kegg data



data(go.sets.hs) #not 100% sure on the difference betwen these sets yet but they seem to retrieve the GO sets
data(go.subs.hs)
#Groups genes based on the 3 gene ontology gene sets
#BP = Biological Process MF = molecular function CC = cellular component
gse16873.bp.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$BP], ref = hn, samp = dcis)
gse16873.mf.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$MF], ref = hn, samp = dcis)
gse16873.cc.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$CC], ref = hn, samp = dcis)

#----------------this is now a more thorougly worked through example----------------



#this is a better way of dividing the samples into their different categories. 
#we will use something like this to pick out the groups we want to compare
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =T)
adh=grep('ADH',cn, ignore.case =T)
dcis=grep('DCIS',cn, ignore.case =T)
print(hn)


#again download sets

data(kegg.gs)
data(go.gs)
lapply(kegg.gs[1:3],head)



kg.hsa=kegg.gsets("hsa") #this picks out the human sets

kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx] #no idea but doesn't seem to work without this step

save(kegg.gs, file="kegg.hsa.sigmet.gsets.RData") #saves the human sets as an R object

#kegg.gsets works with 3000 KEGG species,for examples:
data(korg)
head(korg[,1:3])


#this gets the human sets but for GO instead
go.hs=go.gsets(species="human")
go.bp=go.hs$go.sets[go.hs$go.subs$BP]
go.mf=go.hs$go.sets[go.hs$go.subs$MF]
go.cc=go.hs$go.sets[go.hs$go.subs$CC]
save(go.bp, go.mf, go.cc, file="go.hs.gsets.RData")
#for Bioconductor species supported by go.gsets function:
data(bods)
print(bods)


#These lines are doing the gage and the functions are structured as follows

#name           <- gage(data    ,  genesets used, control group, experimental group)
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)

gse16873.go.p <- gage(gse16873, gsets = go.gs, ref = hn, samp = dcis)

gse16873.bp.p <- gage(gse16873, gsets = go.bp,ref = hn, samp = dcis)



#here we can look at the data
str(gse16873.kegg.p, strict.width='wrap')

head(gse16873.kegg.p$greater[, 1:5], 4)

head(gse16873.kegg.p$less[, 1:5], 4)

head(gse16873.kegg.p$stats[, 1:5], 4)





# this does a two-way analysis
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis, same.dir = F)




write.table(gse16873.kegg.2d.p$greater, file = "gse16873.kegg.2d.p.txt", sep = "\t")

write.table(rbind(gse16873.kegg.p$greater, gse16873.kegg.p$less), file = "gse16873.kegg.p.txt", sep = "\t")



gse16873.kegg.sig<-sigGeneSet(gse16873.kegg.p, outname="gse16873.kegg")

gse16873.kegg.2d.sig<-sigGeneSet(gse16873.kegg.2d.p, outname="gse16873.kegg")



write.table(gse16873.kegg.2d.sig$greater, file = "gse16873.kegg.2d.sig.txt", sep = "\t")


write.table(rbind(gse16873.kegg.sig$greater, gse16873.kegg.sig$less), file = "gse16873.kegg.sig.txt", sep = "\t")





gse16873.kegg.esg.up <- esset.grp(gse16873.kegg.p$greater, gse16873, gsets = kegg.gs, ref = hn, samp = dcis,test4up = T, output = T, outname = "gse16873.kegg.up", make.plot = F)

gse16873.kegg.esg.dn <- esset.grp(gse16873.kegg.p$less, gse16873, gsets = kegg.gs, ref = hn, samp = dcis, test4up = F, output = T, outname = "gse16873.kegg.dn", make.plot = F)


names(gse16873.kegg.esg.up)





#------visualisation----------------------------------------------------------



#heatmap
gs=unique(unlist(kegg.gs[rownames(gse16873.kegg.p$greater)[1:3]]))

essData=essGene(gs, gse16873, ref =hn, samp =dcis)

ref1=1:6

samp1=7:12



for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
   outname = gsub(" |:|/", "_", substr(gs, 10, 100))
   geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
             samp = samp1, outname = outname, txt = T, heatmap = T,
             Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
   }


for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
   outname = gsub(" |:|/", "_", substr(gs, 10, 100))
   outname = paste(outname, "all", sep=".")
   geneData(genes = kegg.gs[[gs]], exprs = gse16873, ref = hn,
              samp = dcis, outname = outname, txt = T, heatmap = T,
              Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}



# interaction networks

gse16873.d <- gse16873[ ,dcis] - gse16873[ ,hn]
path.ids=c("hsa04110 Cell cycle", "hsa00020 Citrate cycle (TCA cycle)")
path.ids2 <- substr(path.ids, 1, 8)
#native KEGG view
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[,1:2], pathway.id = pid, species = "hsa"))








sel <- gse16873.kegg.p$greater[, "q.val"] < 0.1 & !is.na(gse16873.kegg.p$greater[, "q.val"])
path.ids <- rownames(gse16873.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = gse16873.d[,1:2], pathway.id = pid, species = "hsa"))












help(sigGeneSet)
help(gage)
help(pathview)
help(save)