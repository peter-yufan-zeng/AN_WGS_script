#!/usr/bin/env Rscript --vanilla
###This is a helper script to run Chris Mcginis's multise demultiplex script https://github.com/chris-mcginnis-ucsf/MULTI-seq
### run like this Rscript --vanilla /scratch/n/nicholsa/zyfniu/AN_WGS/multiseq_scripts/demultiplex_helper.r /path/read_1 /path/to/read_2
### /path/to/barcodes.Rds /path/to/cell.id.vec.Rds /path/to/save/prefix

library(deMULTIplex)
args = commandArgs(trailingOnly=TRUE)

# test if there are four arguments: if not, return an error
if (length(args)!=5) {
  stop("Five argument must be supplied in this order: /scratch/n/nicholsa/zyfniu/AN_WGS/multiseq_scripts/demultiplex_helper.r /path/read_1 /path/to/read_2
  ### /path/to/barcodes.Rds /path/to/cell.id.vec.Rds /path/to/save/prefix", call.=FALSE)
}
bar.ref <- readRDS(args[3])
cell.id.vec <- readRDS(args[4])
readTable <- MULTIseq.preProcess(R1 = args[1],
                                 R2 = args[2],
                                 cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8)
                                )
saveRDS(readTable,paste0(args[5],"_readTable.rds"))
bar.table <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
saveRDS(bar.table,paste0(args[5],"_barTable.rds"))
