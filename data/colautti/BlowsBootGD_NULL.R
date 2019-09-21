setwd("F://Colautti 2010/2010 Evolution G vs D")
library(boot)

# TO LOAD BOOTSTRAP DATA from previous analysis: 
# BootConst<-read.table("BootConst.txt",header=FALSE)
# BootAngle<-read.table("BootAngle.txt",header=FALSE)
# BootConst<-BootConst[,1]
# BootAngle<-BootAngle[,1]

### IMPORT DATA ###
Fams <- read.table("FamMeans.txt",header=TRUE,row.names=1)
Fams<-Fams[,1:4]
Fams$Lat<-as.factor(Fams$Lat) 
Fams$Fam<-as.factor(Fams$Fam)
Pops <- read.table("PopMeans.txt",header=TRUE,row.names=1) 
Pops<-Pops[,1:3]
Pops$Lat<-as.factor(Pops$Lat) 

### CONVERT LINEAR TRAITS TO SEPARATE COLUMNS ###
FamMeans<-c(1:I(339*14)*0)
FamMeans<-matrix(FamMeans,nrow=339)
FamMeans<-as.data.frame(FamMeans)
names(FamMeans)<-c("Lat","Fam","FDays","FInf","FStemWi","FVeg","HInf","HVeg","Height2","Height4","LnInfW","LnVegW","THeight","TLeafAr")
FamMeans$Lat<-Fams$Lat[Fams$Trait=="FDays"] #note: FDays is arbitrary trait choice here
FamMeans$Fam<-Fams$Fam[Fams$Trait=="FDays"]
for (i in 3:14){
 FamMeans[,i]<-Fams$Estimate[Fams$Trait==names(FamMeans)[i]]
}

PopMeans<-c(1:I(20*13)*0)
PopMeans<-matrix(PopMeans,nrow=20)
PopMeans<-as.data.frame(PopMeans)
names(PopMeans)<-c("Lat","FDays","FInf","FStemWi","FVeg","HInf","HVeg","Height2","Height4","LnInfW","LnVegW","THeight","TLeafAr")
PopMeans$Lat<-Pops$Lat[Pops$Trait=="FDays"] 
for (i in 2:13){
 PopMeans[,i]<-Pops$Estimate[Pops$Trait==names(PopMeans)[i]]
}

### SETUP PARAMETERS FOR BOOTSTRAP ####
nboots<-10000 # number of bootstrap iterations
r<-360/(2*pi)
BootAngle<-0
BootConst<-0
FamVec<-0
FamTraits<-FamMeans[,3:14]
PopTraits<-PopMeans[,2:13]

### Initial estimates of G and D from Mixed Model Data ###
 # PCs of G
FamModel<-prcomp(FamTraits,scale=FALSE)
FamPCs<-FamModel[[2]]
FamPCs<-FamPCs[,1:6]
A<-as.matrix(FamPCs)
At<-t(A)
 # PCs of D
PopModel<-prcomp(PopTraits,scale=FALSE)
PopPCs<-PopModel[[2]]
PopPCs<-PopPCs[,1:6]
B<-as.matrix(PopPCs)
Bt<-t(B)
AtB<-At%*%B
BtA<-Bt%*%A
S<-AtB%*%BtA
PCsS<-as.matrix(eigen(S,symmetric=TRUE))[[2]]
PCsE<-as.matrix(eigen(S,symmetric=TRUE))[[1]]
a1<-as.vector(PCsS[,1])
a2<-as.vector(PCsS[,2])
a3<-as.vector(PCsS[,3])
a4<-as.vector(PCsS[,4])
a5<-as.vector(PCsS[,5])
a6<-as.vector(PCsS[,6])
b1<-A%*%a1
b2<-A%*%a2
b3<-A%*%a3
b4<-A%*%a4
b5<-A%*%a5
b6<-A%*%a6
gmax<-FamPCs[,1]
dmax<-PopPCs[,1]
# Constraint estimated by Krzanowski method
tav<-sum(PCsE)
# Constraint estimated by Schluter method
tang<-acos(abs(gmax%*%dmax))*r


#### NEW BOOTSTRAP Protocol:
for(iter in 1:nboots) {
 BootFam<-FamMeans
 count<-1
 for (i in levels(FamMeans$Lat)){
  temp<-FamMeans[FamMeans$Lat==i,]
  FamNum<-c(1:nrow(temp))
  FamVec<-sample(FamNum,17,replace=TRUE)
   for(j in 1:17){ 
    BootFam[count,]<-temp[FamVec[j],]
    count<-count+1
   }
 }
 ### Calculates Pop Means Based on BLUPs ####
 BootPop<-aggregate(BootFam[,3:14],list(BootFam$Lat),mean)
 names(BootPop)[1]<-"Lat"
 temp<-PopMeans[,2:13]+BootPop[,2:13]
 BootPop<-temp
 temp<-BootFam[,3:14]
 BootFam<-temp
 
 ### Compares G and D and records theta and lambda ####
 FamModel<-prcomp(BootFam,scale=FALSE)
 FamPCs<-FamModel[[2]]
 FamPCs<-FamPCs[,1:6]
 A<-as.matrix(FamPCs)
 At<-t(A) 

 PopModel<-prcomp(BootPop,scale=FALSE)
 PopPCs<-PopModel[[2]]
 PopPCs<-PopPCs[,1:6] 
 B<-as.matrix(PopPCs)
 Bt<-t(B)

 AtB<-At%*%B
 BtA<-Bt%*%A
 S<-AtB%*%BtA
 PCsS<-as.matrix(eigen(S,symmetric=TRUE))[[2]]
 PCsE<-as.matrix(eigen(S,symmetric=TRUE))[[1]]
 a1<-as.vector(PCsS[,1])
 a2<-as.vector(PCsS[,2])
 a3<-as.vector(PCsS[,3])
 a4<-as.vector(PCsS[,4])
 a5<-as.vector(PCsS[,5])
 a6<-as.vector(PCsS[,6])
 b1<-A%*%a1
 b2<-A%*%a2
 b3<-A%*%a3
 b4<-A%*%a4
 b5<-A%*%a5
 b6<-A%*%a6
 gmax<-FamPCs[,1]
 dmax<-PopPCs[,1]

 BootConst[iter]<-sum(PCsE)
 BootAngle[iter]<-acos(abs(gmax%*%dmax))*r

}

# tiff("LambdaBoot.tif",units="in",h=8,w=6,res=200)
# pdf ("LambdaBoot.pdf",width=7.086, height=9.055)
# hist(BootConst,breaks=15)
# hist(BootAngle,breaks=30)



par(oma=c(2,5,2,2),mfrow=c(2,1),mfg=c(1,1),mar=c(4,0,2,0))

hist(BootAngle,main="",xlab="",ylab="",axes=F,breaks=30,xlim=c(90,0))
axis(1,pos=0,at=c(90,60,30,0),las=1,labels=c("90","60","30","0"))
axis(4,pos=0,at=c(0,1600),labels=c("",""),lwd.ticks=0)
axis(3,at=c(0,90),labels=c("",""),lwd.ticks=0)
axis(2,pos=90,at=c(0,400,800,1200,1600,2000,2400),las=1,labels=c("0","400","800","1200","1600","2000",""))

ang<-mean(BootAngle)
#arrows(ang,1200,ang,0,lwd=2,angle=25)
#arrows(tang,1200,tang,0,lwd=2,angle=25,col=rgb(.3,.3,.3),lty=5)

hist(BootConst,main="",xlab="",ylab="",axes=F,breaks=7,xlim=c(0,6))
axis(1,pos=0,at=c(0,2,4,6),las=1,labels=c("0","2","4","6"))
axis(2,pos=0,at=c(0,2000,4000,6000,8000),las=1,labels=c("0","2000","4000","6000","8000"))
axis(3,at=c(0,6),labels=c("",""),lwd.ticks=0)
axis(4,pos=6,at=c(0,8000),labels=c("",""),lwd.ticks=0)
av<-mean(BootConst)
#arrows(av,2000,av,0,lwd=2,angle=25)
#arrows(tav,2000,tav,0,lwd=2,angle=25,col=rgb(.3,.3,.3),lty=5)

par(oma=c(0,0,0,0),mfrow=c(1,1),mfg=c(1,1),mar=c(0,0,0,0))
plot(c(-2,-2),xlim=c(-1,1),ylim=c(-1,1),main="",xlab="",ylab="",axes=F)
text(-0.8,0.95,"a)",cex=2,srt=0)
text(-0.8,-0.05,"b)",cex=2,srt=0)
text(-1,0.1,"Number of bootstrap iterations",cex=2,srt=90)
text(0.1,0.1,expression(theta),font=3,cex=2)
text(0.05,-0.9,substitute(sum()),cex=2)
text(0.13,-0.9,expression(lambda),cex=2)
text(0.17,-0.93,"i",cex=1)

text(-0.7,-1,"No Similarity",cex=1.5)
text(0.9,-1,"Identical",cex=1.5)


# dev.off()

# 95% CI 
low<-(nboots/100)*2.5
high<-(nboots/100)*97.5
BAngle<-sort(BootAngle)
BConst<-sort(BootConst)
angleCI<-c(BAngle[[low]],BAngle[[high]])
indexCI<-c(BConst[[low]],BConst[[high]])

# TO SAVE BOOTSTRAP DATA: write(BootConst,"BootConst.txt",1)
# write(BootAngle,"BootAngle.txt",1)


