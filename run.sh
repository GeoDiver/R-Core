
Rscript download_GEO.R --accession GDS5093 --outrdata ~/Desktop/GDS5093.RData

Rscript dgea.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.RData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --analyse "Boxplot,Volcano,PCA,Heatmap,Clustering" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005 --distance "euclidean" --clustering "average" --heatmaprows 100 --dendrow TRUE --dendcol TRUE --adjmethod fdr --dev TRUE
Rscript dgea_expression.R  --rundir ~/Desktop/ --geneid LOC100288410

Rscript gage.R --accession GDS5093 --dbrdata ~/Desktop/GDS5093.RData --rundir ~/Desktop/ --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control"  --comparisontype ExpVsCtrl --genesettype KEGG --geotype BP --dev TRUE
