setwd("~/BOSTON/Davies/Range expansion/Code/Expected_He/")
#glf <- read.table(file = 'third.het.beagle.gz', header=TRUE)[,-c(1:3)]
glf <- read.table(file = 'lineage_ref.beagle.gz', header=TRUE)[,-c(1:3)]
glf <- array(c(t(as.matrix(glf))), c(3, ncol(glf)/3, nrow(glf)))
EMstep <- function (sfs, GL) rowMeans(prop.table(sfs * GL, 2))

SFS <- matrix(1/3,3,dim(glf)[2])

maxiter <- 200
tol <- 1e-8

for(sample in 1:dim(glf)[2])
{
  for (iter in 1:maxiter)
  {
    upd <- EMstep (SFS[,sample], glf[,sample,])
    if (sqrt(sum((upd - SFS[,sample])^2)) < tol)
      break;
    SFS[,sample] <- upd
  }
  if (iter == maxiter) warning("increase maximum number of iterations")
}



#lineages
print(c('North pop',round(summary(SFS[2,1:94]),4)),quote=F)
print(c('pop2',round(summary(SFS[2,95:115]),4)),quote=F)
print(c('pop4',round(summary(SFS[2,116:158]),4)),quote=F)
#
sfs=SFS[2,1:94]
Site=replicate(94, "Northpop")
Northpop=data.frame(sfs, Site)

sfs=SFS[2,95:115]
Site=replicate(21, "pop2")
pop2=data.frame(sfs, Site)

Site=replicate(43, "pop4")
sfs=SFS[2,116:158]
pop4=data.frame(sfs, Site)

All=rbind(Northpop,pop2,pop4)

library(ggplot2)
# Basic violin plot
colors=c("#ff0088","#caff70", "#006400")
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
p <- ggplot(All, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colors)+
  theme_bw(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 45, vjust = , hjust=1))
p+stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")

t.test(pop2$sfs,pop4$sfs)
#Northpop A
#pop2 B
#pop3 C


#north
print(c('Amakusa',round(summary(SFS[2,1:18]),4)),quote=F)
print(c('Kochi',round(summary(SFS[2,19:37]),4)),quote=F)
print(c('Kushima',round(summary(SFS[2,38:57]),4)),quote=F)
print(c('Kushimoto',round(summary(SFS[2,58:74]),4)),quote=F)
print(c('Shirahama',round(summary(SFS[2,75:92]),4)),quote=F)



sfs=SFS[2,1:18]
Site=replicate(18, "Amakusa")
Amakusa=data.frame(sfs, Site)
sfs=SFS[2,19:37]
Site=replicate(19, "Kochi")
Kochi=data.frame(sfs, Site)
Site=replicate(20, "Kushima")
sfs=SFS[2,38:57]
Kushima=data.frame(sfs, Site)
sfs=SFS[2,58:74]
Site=replicate(17, "Kushimoto")
Kushimoto=data.frame(sfs, Site)
sfs=SFS[2,75:92]
Site=replicate(18, "Shirahama")
Shirahama=data.frame(sfs,Site)
All=rbind(Amakusa,Kochi,Kushima,Kushimoto,Shirahama)

library(ggplot2)
# Basic violin plot
colors=c("#7fabd3ff","#00bfffff", "#6af2ffff","#325fa2","#5087c1")
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
p <- ggplot(All, aes(x=Site, y=sfs,fill=Site)) + 
  geom_violin()+
  scale_fill_manual(values=colors)+
  theme_bw(base_size = 22)+
  ylab(label = "Expected Het")+xlab(label="")+
  theme(axis.text.x = element_text(angle = 45, vjust = , hjust=1))
p+stat_summary(fun.y=mean,size=2, 
               geom="point", color="black")

t.test(Kushimoto$sfs,Kushima$sfs)
Amakusa(AC)
Kochi(B)
Kushima(AB)
Kushimoto(ACD)
Shirahama(D)
#Amakusa 1:18
#Kochi 19:37
#Kushima 38:57
#Kushimoto 58:74
#

###
#Thetas

AmakusaTheta =read.table(file = 'Amakusa.nologs.thetas.idx.raw', header=FALSE)
AmakusaTheta$Site=replicate(nrow(AmakusaTheta), "Amakusa")
KochiTheta=read.table(file = 'Kochi.nologs.thetas.idx.raw', header=FALSE)
KochiTheta$Site=replicate(nrow(KochiTheta), "Kochi")
KushimaTheta=read.table(file = 'Kushima.nologs.thetas.idx.raw', header=FALSE)
KushimaTheta$Site=replicate(nrow(KushimaTheta), "Kushima")
KushimotoTheta=read.table(file = 'Kushimoto.nologs.thetas.idx.raw', header=FALSE)
KushimotoTheta$Site=replicate(nrow(KushimotoTheta), "Kushimoto")
ShirahamaTheta=read.table(file = 'Shirahama.nologs.thetas.idx.raw', header=FALSE)
ShirahamaTheta$Site=replicate(nrow(ShirahamaTheta), "Shirahama")

All=rbind(AmakusaTheta,KochiTheta,KushimaTheta,KushimotoTheta,ShirahamaTheta)

p <- ggplot(All, aes(x=Site, y=V1)) + 
  geom_boxplot()+
  ylab(label = "Expected Het") +ylim(0,0.001)
p

p <- ggplot(All, aes(factor(Site), V1))
p + geom_boxplot(outlier.shape = NA, coef = 0) +ylim(0.0001,0.001) 


