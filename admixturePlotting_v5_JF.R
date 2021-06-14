
setwd("C:/Users/james/Documents/BOSTON/Davies/Range expansion/Code/Admixture/")

# assembling the input table
dir=("C:/Users/james/Documents/BOSTON/Davies/Range expansion/Code/Admixture/") # path to input files
#ngsAdmix
inName="mydata.noLD_k4.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run
#^ this is what I used
#^ does this had no LD? Yes confirmed 12.15.20, moves around a little but same plot. 
#ADMIXTURE
#with LD loci
#inName="myresult2.6.Q"
#inName="myresult.5.Q"





#pops="Adapter2Site.csv" # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.

samples=data.frame(read.csv("C:/Users/james/Documents/BOSTON/Davies/Range expansion/Code/Adapter2Site.csv"))

i2p=samples[-c(2)]
names(i2p)=c("ind","pop")
row.names(i2p)=i2p$ind


bams=data.frame(read.csv("all.bams", header = FALSE)) # list of bam files
colnames(bams)<- "Adapter"
list=bams$Adapter
adapter=gsub(".*Pool","Pool",list)
#adapterall=gsub(".*SRR","SRR",adapter)
#------------

npops=as.numeric(sub(".+(\\d+)\\..+","\\1",inName))

tbl=read.table(paste(dir,inName,sep=""),header=F)
row.names(tbl)= adapter

#new<-merge(tbl,i2p,by="row.names")
new1<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL)
#write.csv(new1, file="bamscl_adaptors.csv")
tbl=new1
head(tbl)
# putting populaitons in desired order (edit pop names as needed or skip to plot them alphabetically)
# tbl$pop=factor(tbl$pop,levels=c("O","K"))
#colors=c("#7fabd3ff","#00bfffff", "#6af2ffff","#325fa2ff","#00ff00ff","#006400FF","#325fa2ff")
colors=c("#ff0085ff","#caff70ff","#660080ff","#006400FF")
#3,2,4,1
tbl$pop=factor(tbl$pop, levels=c("Sekisei","Oura Bay","Kushima","Kochi","Shirahama","Kushimoto","Amakusa"))
#Download below from https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
source("C:/Users/james/Documents/BOSTON/Davies/Range expansion/Code/plot_admixture_v5_function.R")
ords=plotAdmixture(data=tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)



# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.5),collapse=".")) })
cluster.assign=as.data.frame(cluster.admix)
write.csv(cluster.assign, file="cluster.assign.csv")
#write.csv(cluster.assign, file="cluster.assign.HALF.ngsADMIX.csv")
save(cluster.admix,file=paste(inName,"_clusters.RData",sep=""))
