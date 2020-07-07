#!/sw/apps/R/x86_64/3.2.3/milou/bin/Rscript


# Description:
# R-script for running ASCAT

source("AN_WGS/ascat.R")

args = commandArgs(trailingOnly=TRUE)

##args is now a list of character vectors
## First check to see if arguments are passed.
if(length(args)==0){
    stop("No input files supplied\n\nUsage:\nRscript run_ascat.r tumor_baf tumor_logr normal_baf normal_logr\n\n")
} else{
    tumorbaf = args[1]
    tumorlogr = args[2]
    normalbaf = args[3]
    normallogr = args[4]
}

#Load the  data
ascat.bc <- ascat.loadData(Tumor_LogR_file=tumorlogr, Tumor_BAF_file=tumorbaf, Germline_LogR_file=normallogr, Germline_BAF_file=normalbaf)

#Plot the raw data
ascat.plotRawData(ascat.bc)

#Segment the data
ascat.bc <- ascat.aspcf(ascat.bc)

#Plot the segmented data
ascat.plotSegmentedData(ascat.bc)

#Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
ascat.output <- ascat.runAscat(ascat.bc)
