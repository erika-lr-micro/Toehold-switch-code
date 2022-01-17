#Set working directory
setwd("C:/Users/Carlos Z/Downloads")
#Data
toeP = read.table("Tale20.txt", header = FALSE)
str(toeP)
#Sequence from Nupack
seq = read.delim("Tale20_seq.txt", header = FALSE)
#Probability table from Nupack
select = toeP$V2==-1
toeL = toeP[select,]
#Loop to fix table
for(i in 1:3710)
{
  N1 = toeL$V1[i]
  N2 = toeL$V1[i+1]
  if(N2-N1 != 1)
  {
    NN = N1+1
    C2 = -1
    C3 = 0.001
    AnnexData = data.frame("V1"=NN,"V2"=C2,"V3"=C3)
    toeL = rbind(toeL,AnnexData)
    toeL = toeL[order(toeL$V1),]
  }
}
#Output table
Output = data.frame(ID = character(),Seq = character(), P_rna = numeric())
#Loop to select window and calculate total probability
rm(i)
for (i in 1:(length(toeL$V1)-29)) #Loop for segments of 30 bp
{
  Vec = toeL$V3[i:(i+29)]
  Prob = sum(Vec)/30
  ID = paste("T30_",i,"-",i+29, sep="")
  SeqT = substr(as.character(seq$V1[1]),i,i+29)
  AnnexData = data.frame("ID"=ID, "Seq"=SeqT, "P_rna"=Prob)
  Output = rbind(Output,AnnexData)
}

rm(i)
for (i in 1:(length(toeL$V1)-35)) #Loop for segments of 36 bp
{
  Vec = toeL$V3[i:(i+35)]
  Prob = sum(Vec)/36
  ID = paste("T36_",i,"-",i+35, sep="")
  SeqT = substr(as.character(seq$V1[1]),i,i+35)
  AnnexData = data.frame("ID"=ID, "Seq"=SeqT, "P_rna"=Prob)
  Output = rbind(Output,AnnexData)
}

#Function RveComp
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"U"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="U") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}

#Create toehold sensors
Backbone = data.frame(P1=rep("GGG",1), P2=rep("GUUAUAGUUAUGAGACAAGAACAGAGGAGACAUAACAUGAAC",1),P3=rep("GUUAACCUGGCGGCAGCGCAAAAGAUGGGGGAU",1), P4=rep("AAC",1), P5=rep("AACCUGGCGGCAGCGCAAAAG"), P6=rep("AUGGGGGAU",1))
Output$SeqRC = lapply(Output$Seq, rev.comp)
Output$AX = Output$SeqRC
for (i in 1:length(Output$AX))
{
  Value = substr(Output$AX[i], nchar(Output$AX[i])-2, nchar(Output$AX[i]))
  Output$AX[i]=paste(Value,"GUU",sep="")
}
Output$Toehold = paste(Backbone$P1,Output$SeqRC,Backbone$P2,substr(Output$Seq,1,6),Backbone$P4,Output$AX, Backbone$P5,Backbone$P6, sep="")
#Ordering based on the ID so we can cbind the results and be sure they'll be in the same order
Output=Output[order(Output$ID),]
#Packages to download toeholds in different FASTA files
install.packages("dplyr")
library("dplyr")
install.packages("tidyverse")
library(tidyverse)
#Loop to do that multiple times
for (row in 1:nrow(Output)) {
  file.create(paste0(Output[row,1], ".in"))
  fileConn<-file(paste0(Output[row,1], ".in"))
  secuencia = Output[row, 6]
  writeLines(c("1",secuencia,"1"), fileConn)
  close(fileConn)
}
#I have run the files that came out with nupack in the virtualbox with linux
#Files have been saved with the same name plus ".ppairs"
#Creating new data frame
Res_toe = data.frame(ID = character(),P_toe = numeric(), P_sen = numeric())
#Loop for 30nts files
library(ape)
setwd("C:/Users/erika/Escritorio/Universidad/Tesis/Toehold/Archivos_nupack/Resultados30")
files <- list.files(pattern="T*")
for (i in 1:(length(files))){
  filename = substr(files[i], start = 1, stop = nchar(files[i])-7)
  i <- read.table(files[i], header = FALSE, skip=15)
  i <- subset(i, V2==121)
#i refers to the file and j is the variable we use to fix it
  for(j in 1:length(i))
  {
    N1 = i$V1[j]
    N2 = i$V1[j+1]
    if(N2-N1 != 1)
    {
      NN = N1+1
      C2 = 121
      C3 = 0.001
      AnnexData2 = data.frame("V1"=NN,"V2"=C2,"V3"=C3)
      i = rbind(i,AnnexData2)
      i = i[order(i$V1),]
    }
  }
    P_sen=(sum(i$V3))/120
    P_toe=(sum((head(i,33)$V3)))/33
    Res30 = data.frame("ID" = filename,"P_toe" = P_toe, "P_sen" = P_sen)
    Res_toe = rbind(Res_toe,Res30)
  }
#Loop for 36nts files
setwd("C:/Users/erika/Escritorio/Universidad/Tesis/Toehold/Archivos_nupack/Resultados36")
files36 <- list.files(pattern="T*")
for (i in 1:(length(files36))){
  filename36 = substr(files36[i], start = 1, stop = nchar(files36[i])-7)
  i <- read.table(files36[i], header = FALSE, skip=15)
  i <- subset(i, V2==127)
  #i refers to the file and j is the variable we use to fix it
  for(j in 1:length(i))
  {
    N1 = i$V1[j]
    N2 = i$V1[j+1]
    if(N2-N1 != 1)
    {
      NN = N1+1
      C2 = 127
      C3 = 0.001
      AnnexData3 = data.frame("V1"=NN,"V2"=C2,"V3"=C3)
      i = rbind(i,AnnexData3)
      i = i[order(i$V1),]
    }
  }
    P_sen=(sum(i$V3))/126
    P_toe=(sum((head(i,39)$V3)))/39
    Res36 = data.frame("ID" = filename36,"P_toe" = P_toe, "P_sen" = P_sen)
    Res_toe = rbind(Res_toe,Res36)
}
#Organize Res_toe, so it's in the same order as output
Res_toe=Res_toe[order(Res_toe$ID),]
Fdata <- merge(Res_toe, Output)
#Calculating design score
Fdata$Total = (5*Fdata$P_rna) + (4*Fdata$P_toe) + (3*Fdata$P_sen)
#Organizing from highest to lowest
Fdata=Fdata[order(Fdata$Total, decreasing= TRUE),]
#Print table
#Nter = 1-872 
Fdata$Region = ""
for (i in 1:length(Fdata$ID))
{x= Fdata$ID[i]
y = gregexpr(pattern ='-',x)
y = unlist(y)
position = as.numeric(substr(x, y+1, (nchar(x))))
if (position<=871)
{Fdata$Region[i]="Nter"}
else if (position<=2870)
{Fdata$Region[i]="Rep"}
else
{Fdata$Region[i]="Cter"}
}
mean(Fdata$Total)
sd(Fdata$Total)
max(Fdata$Total)
min(Fdata$Total)
DCter<-subset(Fdata, Region=="Cter")
DRep<-subset(Fdata, Region=="Rep")
MaxRep<-subset(Fdata, ID=="T30_2276-2305")
#Setting working directory
setwd("C:/Users/erika/Escritorio/Universidad/Tesis/Toehold/Archivos_nupack")
#Creating a fasta file of the toeholds
file.create("fastarep.FASTA")
for(i in 1:100) {
  row <- DRep[i,]
  cat((paste(">", row$ID, sep="")), file="fastarep.FASTA", fill=TRUE, append=TRUE)
  cat(row$Toehold, file="fastarep.FASTA", fill=TRUE, append=TRUE)
}
