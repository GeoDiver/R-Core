# set up own galaxy

cd /Users/nazrathnawaz/Desktop/R-Core/galaxy
git clone https://github.com/galaxyproject/galaxy/

## copy perl script into Galaxy's tools directory
cd tools
mkdir myTools
cd myTools

## make a tool configuration file - toolExample.xml 
## within the tools/myTools directory


## Make Galaxy aware of the new too
# add these lines to the tool_conf.xml file located in the root directory


## to run 
sh run.sh
# in your browser go to http://localhost:8080







