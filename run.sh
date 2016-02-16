#!/bin/bash 

#### Define Variables ####

ACCESSION='GDS5093'
OUTRDATA=~/Desktop/
FACTOR="disease.state"
DEGUB_MODE=TRUE

POPULATION_A="Dengue Hemorrhagic Fever,Convalescent,Dengue Fever"
POPULATION_B="healthy control"
POPULATION_A_NAME="Dengue"
POPULATION_B_NAME="Normal" 
ANALYSIS_LIST="Boxplot,Volcano,PCA"
DISTANCE_METHOD="euclidean"
CLUSTERING_METHOD="average"
TOPGENE_NO=250
FOLD_CHANGE=0.3
THRESHOLD_VALUE=0.005 
HEATMAP_ROWS=100
DENDOGRAM_ROWS=TRUE
DENDOGRAM_COLUMNS=TRUE
ADJUSTMENT_METHOD=fdr

GENE_ID=LOC100288410

CLUSTER_BY="Complete"
COMPARISON_TYPE=ExpVsCtrl
GENE_SET_TYPE=KEGG
GEO_TYPE=BP

PATH_ID="hsa00480"


#### Run scripts ####

## Install all dependancies
# Rscript Installations.R

## Download and store R data script
Rscript download_GEO.R --accession $ACCESSION --outrdata $OUTRDATA$ACCESSION.RData

## Run overview.R script
Rscript overview.R --accession $ACCESSION --dbrdata $OUTRDATA$ACCESSION.RData --rundir $OUTRDATA/dgea/ --factor $FACTOR --popA $POPULATION_A --popB $POPULATION_B --popname1 $POPULATION_A_NAME --popname2 $POPULATION_B_NAME --analyse $ANALYSIS_LIST --distance $DISTANCE_METHOD --clustering $DISTANCE_METHOD --dev $DEGUB_MODE

## Run DGEA script
Rscript dgea.R --accession $ACCESSION --dbrdata $OUTRDATA$ACCESSION.RData --rundir $OUTRDATA --factor $FACTOR --popA $POPULATION_A --popB $POPULATION_B --popname1 $POPULATION_A_NAME --popname2 $POPULATION_B_NAME --analyse $ANALYSIS_LIST --topgenecount $TOPGENE_NO --foldchange $FOLD_CHANGE --thresholdvalue $THRESHOLD_VALUE --distance $DISTANCE_METHOD --clustering $DISTANCE_METHOD --heatmaprows $HEATMAP_ROWS --dendrow $DENDOGRAM_ROWS --dendcol $DENDOGRAM_COLUMNS --adjmethod $ADJUSTMENT_METHOD --dev $DEGUB_MODE

## Run DGEA Expression level scrip
Rscript dgea_expression.R  --rundir $OUTRDATA --geneid $GENE_ID

## Run GAGE script
Rscript gage.R --accession $ACCESSION --dbrdata $OUTRDATA$ACCESSION.RData --rundir $OUTRDATA --factor $FACTOR --popA $POPULATION_A  --popB $POPULATION_B --comparisontype $COMPARISON_TYPE --genesettype $GENE_SET_TYPE --geotype $GEO_TYPE --distance $DISTANCE_METHOD --clustering $DISTANCE_METHOD --clusterby $CLUSTER_BY --heatmaprows $HEATMAP_ROWS --dendrow $DENDOGRAM_ROWS --dendcol $DENDOGRAM_COLUMNS --dev $DEGUB_MODE

## Run GAGE interaction networks script
Rscript gage_interaction_networks.R --dbrdata $OUTRDATA/kegg.RData --rundir $OUTRDATA --pathid $PATH_ID



