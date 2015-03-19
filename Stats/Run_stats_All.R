
All_L2<-read.delim(file="All_L2", header=TRUE)
All_L3<-read.delim(file="All_L3", header=TRUE)
All_L4<-read.delim(file="All_L4", header=TRUE)
All_L1<-read.delim(file="All_L1", header=TRUE)

#Calculate the aov for ANOVA/TUKEY
AS1<-aov(Percent ~ Sample, data=All_L1)
AS2<-aov(Percent ~ Sample, data=All_L2)
AS3<-aov(Percent ~ Sample, data=All_L3)
AS4<-aov(Percent ~ Sample, data=All_L4)

#Calculate TUKEY HSD
HS1<-HSD.test(AS1, "Sample", group=TRUE)
HS2<-HSD.test(AS2, "Sample", group=TRUE)
HS3<-HSD.test(AS3, "Sample", group=TRUE)
HS4<-HSD.test(AS4, "Sample", group=TRUE)