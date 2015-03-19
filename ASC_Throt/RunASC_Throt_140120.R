RawDataCounts<-read.table(file="SD-NB-Throt.HTSeqCounts.tab", header=TRUE)
source("~/bin/ASC/asc_0.1.6.R")
#Read in the file that should be looped through (Get the Header names for use in the naming of the variables)

DS_out<-"Throt_DS_data_out.rda"
PP_out<-"Throt_PP_data_out.rda"
ci=0
cj=0
#Count variables to determine where we are in the loop
Out=list()

#Loop through the dataframe containing the data just read in.  Assumes that first collumn
#is the protein ids/gene ids for the organisms. 
#Then does a series of comparisons -- where every collumn is compared in succession to the
#columns that follow it. 
label<-as.vector(as.matrix(RawDataCounts[1]))
for (i in RawDataCounts){
  ci=ci+1
  cj=0
  if(ci>1){
    names(i)<-label
    #Adds name vlaues to each of the rows for later comparisons sake
    for (j in RawDataCounts){
      cj=cj+1
      if(cj>ci){
        names(j)<-label
        #Assigns a fabricated variable name (DSDataName) for the DGE approximation and the Posterior Probability estimations for both
        assign(paste("DS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]), sep=""), DGE(X1=i, X2=j))
        #Posterior probability calculation for the delta value off of the raw count data 
        assign(paste("PPS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]),"_2",sep=""), PostProb(d=get(paste("DS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]), sep=""))$delta,x1=i,x2=j,pars=get(paste("DS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]), sep=""))$pars,d0=2))
        assign(paste("PPS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]),"_125",sep=""), PostProb(d=get(paste("DS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]), sep=""))$delta,x1=i,x2=j,pars=get(paste("DS",colnames(RawDataCounts[ci]),"S",colnames(RawDataCounts[cj]), sep=""))$pars,d0=1.25))
        
      }    
    }
  }
save(list=ls(pattern="DS"), file=DS_out)
save(list=ls(pattern="PPS"), file=PP_out)


lst=ls(pattern="PPS")

for (ii in lst){
  filename <- paste(ii, ".txt", sep="")
  print(filename)
  write.table(get(ii), filename, sep="\t",quote=FALSE)
}
}

