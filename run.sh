#!/bin/bash 

#### Define Variables ####

ACCESSION='GDS5093'
OUTRDATA='analysis/'
FACTOR="disease.state"
DEGUB_MODE=TRUE

POPULATION_A="Dengue Hemorrhagic Fever,Convalescent,Dengue Fever"
POPULATION_B="healthy control"
POPULATION_A_NAME="Dengue"
POPULATION_B_NAME="Normal" 
ANALYSIS_LIST="Boxplot,Volcano,PCA,Heatmap"
DISTANCE_METHOD="euclidean"
CLUSTERING_METHOD="average"
TOPGENE_NO=250
FOLD_CHANGE=0.0
THRESHOLD_VALUE=0.005 
HEATMAP_ROWS=100
DENDOGRAM_ROWS=TRUE
DENDOGRAM_COLUMNS=TRUE
ADJUSTMENT_METHOD=fdr

GENE_ID=LOC100288410

CLUSTER_BY="Complete"
COMPARISON_TYPE="ExpVsCtrl"
GENE_SET_TYPE="KEGG"
GEO_TYPE="BP"

PATH_ID="hsa00480"


#### Run scripts ####
function check_if_file_exists {
  if [ -s $1 ]; then 
    echo "*** File Created: $1"
  else
    echo "*** ERROR: File Does Not exist: $1"
    exit 1
  fi
}

mkdir $OUTRDATA

##############################
#### download_GEO.R
##############################
## Download and store R data script
Rscript download_GEO.R --accession "$ACCESSION" --outrdata $OUTRDATA/$ACCESSION.RData
check_if_file_exists $OUTRDATA/$ACCESSION.RData
echo
echo "###################"
echo "###################"
echo
echo
##############################
#### Overview.R
##############################
echo Rscript overview.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA" --factor $FACTOR --popA "$POPULATION_A" --popB "$POPULATION_B" --popname1 "$POPULATION_A_NAME" --popname2 "$POPULATION_B_NAME" --analyse "$ANALYSIS_LIST" --dev "$DEGUB_MODE"
Rscript overview.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA" --factor $FACTOR --popA "$POPULATION_A" --popB "$POPULATION_B" --popname1 "$POPULATION_A_NAME" --popname2 "$POPULATION_B_NAME" --analyse "$ANALYSIS_LIST" --dev "$DEGUB_MODE"
echo
echo "## Verifying Output"
check_if_file_exists $OUTRDATA/boxplot.png
check_if_file_exists $OUTRDATA/data.json
echo
echo "###################"
echo "###################"
echo
echo

##############################
#### dgea.R
##############################
echo Rscript dgea.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA" --factor "$FACTOR" --popA "$POPULATION_A" --popB "$POPULATION_B" --popname1 "$POPULATION_A_NAME" --popname2 "$POPULATION_B_NAME" --analyse "$ANALYSIS_LIST" --topgenecount "$TOPGENE_NO" --foldchange "$FOLD_CHANGE" --thresholdvalue "$THRESHOLD_VALUE" --distance "$DISTANCE_METHOD" --clustering "$CLUSTERING_METHOD" --clusterby "$CLUSTER_BY" --heatmaprows "$HEATMAP_ROWS" --dendrow "$DENDOGRAM_ROWS" --dendcol "$DENDOGRAM_COLUMNS" --adjmethod "$ADJUSTMENT_METHOD" --dev "$DEGUB_MODE"
Rscript dgea.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA" --factor "$FACTOR" --popA "$POPULATION_A" --popB "$POPULATION_B" --popname1 "$POPULATION_A_NAME" --popname2 "$POPULATION_B_NAME" --analyse "$ANALYSIS_LIST" --topgenecount "$TOPGENE_NO" --foldchange "$FOLD_CHANGE" --thresholdvalue "$THRESHOLD_VALUE" --distance "$DISTANCE_METHOD" --clustering "$CLUSTERING_METHOD" --clusterby "$CLUSTER_BY" --heatmaprows "$HEATMAP_ROWS" --dendrow "$DENDOGRAM_ROWS" --dendcol "$DENDOGRAM_COLUMNS" --adjmethod "$ADJUSTMENT_METHOD" --dev "$DEGUB_MODE"
echo
echo "## Verifying Output"
check_if_file_exists $OUTRDATA/dgea_data.json
check_if_file_exists $OUTRDATA/dgea_heatmap.svg
check_if_file_exists $OUTRDATA/dgea_toptable.RData
check_if_file_exists $OUTRDATA/dgea_toptable.tsv
check_if_file_exists $OUTRDATA/dgea_volcano.png
echo
echo "###################"
echo "###################"
echo
echo

##############################
#### dgea_expression.R
##############################
echo Rscript dgea_expression.R  --rundir "$OUTRDATA/" --geneid "$GENE_ID"
Rscript dgea_expression.R  --rundir "$OUTRDATA/" --geneid "$GENE_ID"
echo
echo "## Verifying Output"
check_if_file_exists $OUTRDATA/dgea_$GENE_ID.json
echo
echo "###################"
echo "###################"
echo
echo

##############################
#### gage.R
##############################
echo Rscript gage.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA/" --factor "$FACTOR" --popA "$POPULATION_A"  --popB "$POPULATION_B" --comparisontype "$COMPARISON_TYPE" --genesettype "$GENE_SET_TYPE" --distance "$DISTANCE_METHOD" --clustering "$DISTANCE_METHOD" --clusterby "$CLUSTER_BY" --heatmaprows "$HEATMAP_ROWS" --dendrow "$DENDOGRAM_ROWS" --dendcol "$DENDOGRAM_COLUMNS" --dev "$DEGUB_MODE"
Rscript gage.R --accession "$ACCESSION" --dbrdata "$OUTRDATA/$ACCESSION.RData" --rundir "$OUTRDATA/" --factor "$FACTOR" --popA "$POPULATION_A"  --popB "$POPULATION_B" --comparisontype "$COMPARISON_TYPE" --genesettype "$GENE_SET_TYPE" --distance "$DISTANCE_METHOD" --clustering "$DISTANCE_METHOD" --clusterby "$CLUSTER_BY" --heatmaprows "$HEATMAP_ROWS" --dendrow "$DENDOGRAM_ROWS" --dendcol "$DENDOGRAM_COLUMNS" --dev "$DEGUB_MODE"
echo
echo "## Verifying Output"
check_if_file_exists $OUTRDATA/gage_data.json
check_if_file_exists $OUTRDATA/gage.RData
check_if_file_exists $OUTRDATA/gage_heatmap.svg
check_if_file_exists $OUTRDATA/gage_toptable.tsv
echo
echo "###################"
echo "###################"
echo
echo

##############################
#### gage_interaction_networks.R
##############################
cd $OUTRDATA
echo Rscript ../gage_interaction_networks.R --rundir ./ --pathid $PATH_ID
Rscript ../gage_interaction_networks.R --rundir ./ --pathid $PATH_ID
echo
echo "## Verifying Output"
check_if_file_exists $PATH_ID.gage_pathway.multi.png
check_if_file_exists $PATH_ID.png
check_if_file_exists $PATH_ID.xml
echo
echo "###################"
echo "###################"
echo
echo

