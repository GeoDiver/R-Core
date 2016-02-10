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
biocLite(c("GEOquery","gage","gageData","GO.db", "pathview"))
install.packages('argparser')
install.packages('Cairo')
install.packages('dendextend')
install.packages('DMwR')
install.packages('GEOquery')
install.packages('ggplot2')
install.packages('gplots')
install.packages('jsonlite')
install.packages('limma')
install.packages('pheatmap')
install.packages('plyr')
install.packages('RColorBrewer')
install.packages('reshape2')
install.packages('squash')