#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : Installations.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Pre-Installation libraries
# ---------------------------------------------------------

#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

source('http://bioconductor.org/biocLite.R')
biocLite("GEOquery")
biocLite("gage")
biocLite("gageData")
biocLite("GO.db")
biocLite("pathview")
biocLite("limma")
biocLite("org.Mm.eg.db")
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
# DMwR          - https://cran.r-project.org/web/packages/DMwR/index.html
# ggplot2       - https://cran.r-project.org/web/packages/ggplot2/index.html
# gplots        - https://cran.r-project.org/web/packages/gplots/index.html
# jsonlite      - https://cran.r-project.org/web/packages/jsonlite/index.html
# pheatmap      - https://cran.r-project.org/web/packages/pheatmap/index.html
# plyr          - https://cran.r-project.org/web/packages/plyr/index.html
# RColorBrewer  - https://cran.r-project.org/web/packages/RColorBrewer/index.html
# reshape2      - https://cran.r-project.org/web/packages/reshape2/index.html
# squash        - https://cran.r-project.org/web/packages/squash/index.html
