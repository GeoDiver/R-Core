#!/usr/bin/Rsript
# ---------------------------------------------------------------------------
# Filename      : DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Retrieve individual gene expressions and convert to JSON
# Run           : Rscript dgea_expression.R  --rundir ~/Desktop/ --geneid LOC100288410
# ---------------------------------------------------------------------------

suppressMessages(library('argparser'))    # Argument passing
suppressMessages(library('jsonlite'))     # Convert R object to JSON format

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--rundir", help = "Full file path of rData file") 
parser <- add_argument(parser, "--geneid", help = "Row Id of the X matrix")    

# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                          Loading Saved Dataset                            #
#############################################################################

filename <- paste(argv$rundir,"dgea_toptable.RData", sep = "")

if (file.exists(filename)){
    load(file = filename)
}else {
    # Exit with error code 1
    print("ERROR:File not found.Run the analysis first to see Expression Profile")
    quit(save = "no", status = 1, runLast = FALSE)
}

#############################################################################
#                          Retrieve Expression data                         #
#############################################################################

# check directory exists and availability of toptable expression data
if ((!is.na(argv$rundir)) && (!is.na(X.toptable))) {
    
    # get indexes of group 1
    index.group1 <- which((expression.info["population"] == "Group1") == TRUE)
    # get sample names (x) and expression data (y)
    g1 <- list( x = names(X.toptable[argv$geneid, index.group1]),
                y = as.double(X.toptable[argv$geneid, index.group1]))
  
    # get indexes of group 2
    index.group2 <- which((expression.info["population"] == "Group2") == TRUE)
    # get sample names (x) and expression data (y)
    g2 <- list(x = names(X.toptable[argv$geneid, index.group2]),
                   y = as.double(X.toptable[argv$geneid, index.group2]))
    
    # write to a file
    filename <- paste(argv$rundir,"dgea_",argv$geneid,".json", sep = "")
    write(toJSON(list(group1 = g1 , group2 = g2)), filename)
} else {
    # Exit with error code 1
    print("ERROR:File not found or toptable data not found!")
    quit(save = "no", status = 1, runLast = FALSE)
}
