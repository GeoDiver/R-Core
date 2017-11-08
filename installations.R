#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : Installations.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Pre-Installation libraries
# ---------------------------------------------------------

#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("gage")
biocLite("gageData")
biocLite("GO.db")
biocLite("pathview")
biocLite("limma")
biocLite("impute")
# Species database
biocLite("org.Ag.eg.db") # "Anopheles" "aga" "eg"
biocLite("org.At.tair.db") # "Arabidopsis" "ath" "tair"
biocLite("org.Bt.eg.db") # "Bovine" "bta" "eg"
biocLite("org.Ce.eg.db") # "Worm" "cel" "eg"
biocLite("org.Cf.eg.db") # "Canine" "cfa" "eg"
biocLite("org.Dm.eg.db") # "Fly" "dme" "eg"
biocLite("org.Dr.eg.db") # "Zebrafish" "dre" "eg"
biocLite("org.EcK12.eg.db") # "E coli strain K12" "eco" "eg"
biocLite("org.EcSakai.eg.db") # "E coli strain Sakai" "ecs" "eg"
biocLite("org.Gg.eg.db") # "Chicken" "gga" "eg"
biocLite("org.Hs.eg.db") # "Human" "hsa" "eg"
biocLite("org.Mm.eg.db") # "Mouse" "mmu" "eg"
biocLite("org.Mmu.eg.db") # "Rhesus" "mcc" "eg"
biocLite("org.Pf.plasmo.db") # "Malaria" "pfa" "orf"
biocLite("org.Pt.eg.db") # "Chimp" "ptr" "eg"
biocLite("org.Rn.eg.db") # "Rat" "rno" "eg"
biocLite("org.Sc.sgd.db") # "Yeast" "sce" "orf"
biocLite("org.Ss.eg.db") # "Pig" "ssc" "eg"
biocLite("org.Xl.eg.db") # "Xenopus" "xla" "eg"
install.packages("argparser")
install.packages("Cairo")
install.packages("dendextend")
install.packages("DMwR")
install.packages("ggplot2")
install.packages("gplots")
install.packages("jsonlite")
install.packages("pheatmap")
install.packages("plyr")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("squash")

#############################################################################
#            For additional information on required packages                #
#############################################################################

# GEOquery      - http://bioconductor.org/packages/release/bioc/html/GEOquery.html
# gage          - http://bioconductor.org/packages/release/bioc/html/gage.html
# gageData      - https://bioconductor.org/packages/release/data/experiment/html/gageData.html
# GO.db         - https://bioconductor.org/packages/release/data/annotation/html/GO.db.html
# pathview      - http://bioconductor.org/packages/release/bioc/html/pathview.html
# limma         - https://bioconductor.org/packages/release/bioc/html/limma.html
# org.Mm.eg.db  - http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
# argparser     - https://cran.r-project.org/web/packages/argparser/index.html
# Cairo         - https://cran.r-project.org/web/packages/Cairo/index.html
# dendextend    - https://cran.r-project.org/web/packages/dendextend/index.html
# impute        - http://svitsrv25.epfl.ch/R-doc/library/impute/html/impute.knn.html
# DMwR          - https://cran.r-project.org/web/packages/DMwR/index.html
# ggplot2       - https://cran.r-project.org/web/packages/ggplot2/index.html
# gplots        - https://cran.r-project.org/web/packages/gplots/index.html
# jsonlite      - https://cran.r-project.org/web/packages/jsonlite/index.html
# pheatmap      - https://cran.r-project.org/web/packages/pheatmap/index.html
# plyr          - https://cran.r-project.org/web/packages/plyr/index.html
# RColorBrewer  - https://cran.r-project.org/web/packages/RColorBrewer/index.html
# reshape2      - https://cran.r-project.org/web/packages/reshape2/index.html
# squash        - https://cran.r-project.org/web/packages/squash/index.html
