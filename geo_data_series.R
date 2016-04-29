
#### GeoDiver 2.0 for GSE Data Series ####


#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

library('Cairo')
library('dendextend')
library('GEOquery')
library('ggplot2')
library('gplots')
library('rjson')
library('limma')
library('plyr')
library('RColorBrewer')
library('reshape2')


#############################################################################
#                        Testing Variables                         #
#############################################################################

output.dir     <- "/Users/nazrathnawaz/Desktop/"

# --------- Geo DataSet Input ------------

accession.id    <- 'GSE29290'
factor.type     <- 'disease.state'
population1     <- c('Dengue Hemorrhagic Fever','Convalescent')
population2     <- c('healthy control')
pop.name1       <- "Dengue"
pop.name2       <- "Normal"
pop.colour1     <- "#b71c1c" # Red
pop.colour2     <- "#0d47a1" # Blue

# --------- Volcano Plot ------------
no.of.top.genes <- 250
toptable.sortby <- "p"
fold.change 	<- 0.3
threshold.value <- 0.005


#############################################################################
#                        GEO Input                          #
#############################################################################

# import data sets and process into expression data

accession.id    	<- 'GDS5093'
gse  				<- getGEO(accession.id, GSEMatrix = TRUE)       # Load GEO data
met              	<- Meta(gse)
eset             	<- GDS2eSet(gse, do.log2=TRUE)                  # Convert into ExpressionSet Object
X                	<- exprs(eset)                                  # Get Expression Data
gene.names        	<- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
names(gene.names) 	<- rownames(X)
pClass           	<- pData(eset)[factor.type]
colnames(pClass) 	<- 'factor.type'
samples          	<- rownames(pClass)



accession.id    	<- 'GSE1871'
gse              	<- getGEO(accession.id, GSEMatrix = TRUE)       # Load GEO data

if (length(gse) > 1) idx <- grep("GPL1261", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]

getGSEDataTables(gse)

parseGDS(gse)

#############################################################################
#                        Two Population Preparation                          #
#############################################################################

# Create a data frame with the factors
expression.info  <- data.frame(pClass, Sample = samples, row.names = samples)

# Introduce two columns to expression.info -
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
expression.info <- within(expression.info, {
    population        = ifelse (factor.type %in% population1,'Group1', # if true
                                ifelse( factor.type %in% population2, 'Group2', NA) ) # if false
    population.colour = ifelse (factor.type %in% population1,pop.colour1, # if true
                                ifelse( factor.type %in% population2, pop.colour2, '#000000') ) # if false
})

# Convert to a factor
expression.info$population <- as.factor(expression.info$population)

data <- within(melt(X), {
    phenotypes = expression.info[Var2, 'factor.type']
})

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info),]
X <- X[,(colnames(X) %in% rownames(expression.info))]

# Created a new Phenotype class
newPClass           <- expression.info$population
names(newPClass)    <- expression.info$Sample

#############################################################################
#                        Meta Data Access                         #
#############################################################################

# Obtain a vector of the possible factors
get.factors <- function(pDat, remove.list = c('sample', 'description', 'individual') ) {
    return ( colnames(pDat)[ !colnames(pDat) %in% remove.list ] )
}

generate.geo.summary.json <- function(pDat, met, factors)  { 
    if (missing(factors)) {
        factors <- get.factors(pDat)
    }


    results.data <- list(Title = met$title,
                         Description = met$description[1],
                         Sample_Count = met$sample_count,
                         Feature_Count = met$feature_count,
                         Sample_Organism = met$sample_organism,
                         Platform = met$platform,
                         Platform_Technology_Type = met$platform_technology_type,
                         Sample_Type = met$sample_type,
                         Factor = lapply(pDat[factors], levels),
                         Pubmed_ID = met$pubmed_id,
                         Reference = met$ref,
                         Reference_Series = met$reference_series,
                         Update_Date = met$update_date)
    data.JSON <- toJSON( results.data )
    return ( data.JSON )
}


#############################################################################
#                        Top Table                          #
#############################################################################

find.toptable <- function(X, newPClass, toptable.sortby, no.of.top.genes, gene.names){

    design  <- model.matrix(~0 + newPClass)

    # plots linear model for each gene and estimate fold changes and standard errors
    fit     <- lmFit(X, design)

    # set contrasts for all classes
    contrasts <- makeContrasts(contrasts="newPClassGroup1-newPClassGroup2",
                               levels = design)

    fit <- contrasts.fit(fit, contrasts)

    # empirical Bayes smoothing to standard errors
    fit <- eBayes(fit)

    # Sort.by shoudl be variable - 'p' or 'LogFC'
    toptable <- topTable(fit, sort.by= toptable.sortby, number=no.of.top.genes, genelist = gene.names)

    return(toptable)
}


filtered.toptable <- function(toptable, gene.names){
    toptable['gene'] <- gene.names[rownames(toptable)]

    #Create Sub Data
    X.toptable <- X[rownames(toptable),]
    rownames(X.toptable) <- gene.names[rownames(toptable)]

    return(X.toptable)
}
#############################################################################
#                        Graphical Representation                          #
#############################################################################

samples.boxplot <- function(){
    boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = phenotypes)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = as.vector(expression.info$population.colour)), legend.position = 'right')+ labs(x = 'Samples', y = 'Expression Levels')
    # store Boxplot as an .svg file in the working directory
    filename <- paste(output.dir,"boxplot.svg",sep = "")
    ggsave(filename, plot=boxplot, width = 8, height = 4)
}

heatmap <- function(X, sample.colours, cv = TRUE, rv = TRUE){
    # store Heatmap as an .svg file in the working directory
    filename <- paste(output.dir,"Heatmap.svg",sep = "")
    CairoSVG(file = filename)
    color_scale <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
    heatmap1 <- heatmap.2(X, col=color_scale, scale='row',
                          key=T, keysize=1,
                          dendrogram='column', density.info='none',
                          trace='none', cexCol=0.5, cexRow=0.1,
                          ColSideColors = sample.colours,
                          Colv = cv, Rowv = rv)
    dev.off()
}

adj.p.val.histogram <- function(toptable){
    # store Histogram as an .svg file in the working directory
    filename <- paste(output.dir,"Histogram.svg",sep = "")
    CairoSVG(file = filename)
    hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', xlab = "Adjusted p-values", main=NULL)
    dev.off()
}

volcanoplot1 <- function(toptable){
    # store Volcano plot as an .svg file in the working directory
    filename <- paste(output.dir,"Volcano.svg",sep = "")
    CairoSVG(file = filename)
    with(toptable, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim = c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1),ylim = c(0, max(-log10(toptable$P.Value))+0.5)))
    #volcanoplot(fit, coef=1, highlight=20, names=gene.names, col='steelblue', xlab='Log Fold Change',
    #            ylab='Log Odds', pch=16, cex=0.5)
    dev.off()
}


                                                #Bonferroni cut-off    
volcanoplot2 <- function(toptable,fold.change, t = 0.05/length(gene.names)){

    # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    toptable$threshold = as.factor(abs(toptable$logFC) > fold.change & toptable$P.Value < t)

    # Construct the plot object
    vol = ggplot(data=toptable, aes(x=toptable$logFC, y=-log10(toptable$P.Value), colour=threshold)) +
        geom_point(alpha=0.4, size=1.75)  + xlim(c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1)) + ylim(c(0, max(-log10(toptable$P.Value))+0.5)) +
        xlab("log2 fold change") + ylab("-log10 p-value")
    filename <- paste(output.dir,"Volcanoplot2.svg",sep = "")
    ggsave(filename, plot=vol, height = 6, width = 6)
}

top.genes <- function(toptable, n){
    # store Top20 genes as a .csv file in the working directory
    filename <- paste(output.dir,"TopGenes.csv",sep = "")
    write.csv(toptable[1:n,], file = filename)
}

clustering <- function(dist.method = "euclidean", clust.method = "average"){
    hc <- hclust(dist(t(X),dist.method), clust.method)
    dend <- as.dendrogram(hc)
    labels_colors(dend) <- expression.info$population.colour[order.dendrogram(dend)]
    filename <- paste(output.dir,"Cluster.svg",sep = "")
    CairoSVG(file = filename)
    plot(dend, main = "Cluster Dendrogram", xlab = "Samples")
    dev.off()
}

# ------- Principal Component Analysis ----------
#
# Xpca <- prcomp(t(X), scale= TRUE)
# attributes(Xpca)
#
# s <- summary(Xpca)
# str(s)
#
# # Explained Variance
# expVar <- s$importance[2,] * 100   # convert to %
# expVar
#
# barplot(expVar, xlab = "Principal Components", ylab = "Explained Variance (%)",
#         col ="steelblue")
#
# # Cumulative Variance
# cumVar <- s$importance[3,] * 100
#
# plot(cumVar, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance",
#      xlab="Principal Components")
#
# Xscores <- Xpca$x
#
# Xloadings <- Xpca$rotation
#
#
# pValues <- c()
# for (i in 1:nrow(X)) {
#     res <- t.test(X[i,controlIndex],X[i,virusIndex])
#     pValues <- c(pValues, res$p.value)
# }
# names(pValues) <- geneNames
# topGenesIndx <- order(pValues)[1:50]
# topGeneNames <- geneNames[topGenesIndx]
# topGeneNames
#
# Xpca <- prcomp(t(X[topGenesIndx, ]), scale= TRUE)
# # Explained Variance
# expVar <- s$importance[2,] * 100
# expVar
#
# barplot(expVar, xlab = "Principal Components", ylab = "Explained Variance (%)",
#         col ="steelblue")
#
# # Cumulative Variance
# cumVar <- s$importance[3,] * 100
# cumVar
#
# plot(cumVar, type="o" ,col="black", pch=21, bg="blue", cex=0.8, ylab="Cumulative Variance",
#      xlab="Principal Components")

# Xscores <- Xpca$x
# plot(Xscores, xlab="PC1", ylab="PC2", pch=21, bg=colour, cex=0.7,
#      cex.lab=0.7, cex.axis = 0.7)

#############################################################################
#                        Function Calling                          #
#############################################################################

exportJson <- generate.geo.summary.json(pData(eset), met)
write(exportJson, 'factors.json')

no.of.top.genes<-250
samples.boxplot()
toptable <- find.toptable(X, newPClass, toptable.sortby, no.of.top.genes, gene.names)
X.toptable <- filtered.toptable(toptable, gene.names)
heatmap(X.toptable,expression.info$population.colour)
adj.p.val.histogram(toptable)
toptable.all <- find.toptable(X, newPClass, toptable.sortby, length(gene.names) , gene.names)
volcanoplot1(toptable)
volcanoplot2(toptable.all,fold.change)
top.genes(toptable,no.of.top.genes)
clustering()


