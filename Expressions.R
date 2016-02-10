#!/usr/bin/Rsript
# ---------------------------------------------------------------------------
# Filename      : DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Retrieve individual gene expressions and convert to JSON
# Run           : Rscript Expressions.R  --rundir ~/Desktop/ --geneid LOC100288410
# ---------------------------------------------------------------------------

library('argparser')    # Argument passing
library('jsonlite')     # Convert R object to JSON format

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--rundir", help="Full file path of rData file") 
parser <- add_argument(parser, "--geneid", help="Row Id of the X matrix")    

# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                          Loading Saved Dataset                            #
#############################################################################

filename <- paste(argv$rundir,"expressionprofile.rData", sep = "")

if (file.exists(filename)){
    load(file = filename)
}else{
    print("ERROR:File not found")
    q(save = "no")
}

#############################################################################
#                          Retrieve Expression data                         #
#############################################################################

if((!is.na(argv$rundir))&&(!is.na(X.toptable))){
    index.group1 <- which((expression.info['population']== 'Group1')==TRUE)
    g1 <- list( x = names(X.toptable[argv$geneid, index.group1]),
                y = as.double(X.toptable[argv$geneid, index.group1]))
  
    index.group2 <- which((expression.info['population']== 'Group2')==TRUE)
    g2 <- list(x = names(X.toptable[argv$geneid, index.group2]),
                   y = as.double(X.toptable[argv$geneid, index.group2]))
    filename <- paste(argv$rundir,argv$geneid,".json", sep = "")
    write(toJSON(list(group1 = g1 ,group2 = g2)), filename)
}