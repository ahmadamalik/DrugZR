# DrugZR


Usage:

CountTab<-read.table("DS000015315_drugZ_countTable.txt",header=TRUE)


treatment<-c("MHCI.hi_NA_1","MHCI.hi_NA_2")
control<-c("Unsorted_NA_1","Unsorted_NA_2")


DrugZRes<-DrugZR(se=countTab,control,treatment,min_observations=1,pseudocount=5,seTrue=FALSE)
