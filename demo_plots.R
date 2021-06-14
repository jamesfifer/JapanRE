setwd("~/BOSTON/Davies/Range expansion/Code/Demographic_Analysis/Final_Runs")

pop24=read.table(file="boot.pop2pop4.out.params")
nepop24=pop24[,8:11]
Tpop24=pop24[,12:14]
m_2=pop24[,1]+pop24[,2]
mig24_3=cbind(pop24[,3:6],m_2)

mig24_3_meld=melt(data=mig24_3, measure.vars=colnames(mig24_3))
plot24=ggplot(mig24_3_meld, aes(x=factor(variable, level=c('m_2','m12_3','m21_3','m12_3i','m21_3i')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("migration")+xlab(NULL)+ theme(axis.text.x = element_text(angle = 45)) 
plot24
library(reshape2)
library(ggplot2)

statpopne24=melt(data=nepop24, measure.vars=colnames(nepop24))

ggplot(statpopne24, aes(x=factor(variable, level=c('nu1_2','nu2_2','nu1_3','nu2_3')), y=log(value)))+
  geom_boxplot()+
  scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)

statpopT24=melt(data=Tpop24, measure.vars=colnames(Tpop24))
ggplot(statpopT24,aes(x=factor(variable, level=c('T0','T1','T2')),y=value))+
geom_boxplot()+
  scale_y_continuous(name="Years from present", labels = scales::comma,limits=c(0,300000))+
  theme_classic(base_size = 22)+xlab(NULL) 
#
pop2N=read.table(file="boot.pop2popN.out.params")
nepop2N=pop2N[,8:11]
Tpop2N=pop2N[,12:14]
m_2=pop2N[,1]+pop2N[,2]
mig24_3=cbind(pop2N[,3:6],m_2)
mig24_3_meld=melt(data=mig24_3, measure.vars=colnames(mig24_3))
plot2N=ggplot(mig24_3_meld, aes(x=factor(variable, level=c('m_2','m12_3','m21_3','m12_3i','m21_3i')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("migration")+xlab(NULL)+ theme(axis.text.x = element_text(angle = 45))
plot2N
library(reshape2)
library(ggplot2)

library(reshape2)
library(ggplot2)

statpopne2N=melt(data=nepop2N, measure.vars=colnames(nepop2N))

ggplot(statpopne2N, aes(x=factor(variable, level=c('nu1_2','nu2_2','nu1_3','nu2_3')), y=log(value)))+
  geom_boxplot()+
  scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)

statpopT2N=melt(data=Tpop2N, measure.vars=colnames(Tpop2N))
ggplot(statpopT2N,aes(x=factor(variable, level=c('T0','T1','T2')),y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Years from present", labels = scales::comma,limits=c(0,300000))+
  theme_classic(base_size = 22)+xlab(NULL) 

#
pop4N=read.table(file="boot.pop4popN.out.params")
nepop4N=pop4N[,5:8]
Tpop4N=pop4N[,9:10]

m_2=pop4N[,1:4]
mig_2meld=melt(data=m_2, measure.vars=colnames(m_2))
plot4N=ggplot(mig_2meld, aes(x=factor(variable, level=c('m12_2','m21_2','m12_2i','m21_2i')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("migration")+xlab(NULL)+ theme(axis.text.x = element_text(angle = 45))
plot4N
library(reshape2)
library(ggplot2)

statpopne4N=melt(data=nepop4N, measure.vars=colnames(nepop4N))

ggplot(statpopne4N, aes(x=factor(variable, level=c('nu1_1','nu2_1','nu1_2','nu2_2')), y=log(value)))+
  geom_boxplot()+
  scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)

statpopT4N=melt(data=Tpop4N, measure.vars=colnames(Tpop4N))
ggplot(statpopT4N,aes(x=factor(variable, level=c('T0','T1')),y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Years from present", labels = scales::comma,limits=c(0,300000))+
  theme_classic(base_size = 22)+xlab(NULL) 
#all migs
library(cowplot)
plot_grid(
  plot24,plot2N,NULL,plot4N,
  ncol = 2,
  label_size = 12,
  align = "v"
)

###########################################################################################################
#NORTH

AmaKusha=read.table(file="boot.AmakusaKushima.out.params")
mean(AmaKusha$m12)
mean(AmaKusha$m21)
median(AmaKusha$m12)
median(AmaKusha$m21)
AmaKushameanmig12=mean(AmaKusha$m12)
AmaKushameanmig21=mean(AmaKusha$m21)
AmaKushaNe1=AmaKusha$nu1_1
AmaKushaNe2=AmaKusha$nu2_1

KocKusha=read.table(file="boot.KochiKushima.out.params")
mean(KocKusha$m)
median(KocKusha$m)
KocKushameanmig=mean(KocKusha$m)
KocKushaNe1=KocKusha$nu1_1
KocKushaNe2=KocKusha$nu2_1

KocKusho=read.table(file="boot.KochiKushimoto.out.params")
mean(KocKusho$m12)
mean(KocKusho$m21)
median(KocKusho$m12)
median(KocKusho$m21)
KocKushomeanmig12=mean(KocKusho$m12)
KocKushomeanmig21=mean(KocKusho$m21)
KocKushoNe1=KocKusho$nu1_1
KocKushoNe2=KocKusho$nu2_1

KocShi=read.table(file="boot.KochiShirahama.out.params")
mean(KocShi$m12)
mean(KocShi$m21)
median(KocShi$m12)
median(KocShi$m21)
KocShimeanmig12=mean(KocShi$m12)
KocShimeanmig21=mean(KocShi$m21)
KocShiNe1=KocShi$nu1_1
KocShiNe2=KocShi$nu2_1

KushaKusho=read.table(file="boot.KushimaKushimoto.out.params")
mean(KushaKusho$m12_2)
mean(KushaKusho$m21_2)
median(KushaKusho$m12_2)
median(KushaKusho$m21_2)
KushaKushomeanmig12=mean(KushaKusho$m12_2)
KushaKushomeanmig21=mean(KushaKusho$m21_2)
KushaKushoNe1=KushaKusho$nu1_2
KushaKushoNe2=KushaKusho$nu2_2

KushaShi=read.table(file="boot.KushimaShirahama.out.params")
mean(KushaShi$m12_2)
mean(KushaShi$m21_2)
median(KushaShi$m12_2)
median(KushaShi$m21_2)
KushaShimeanmig12=mean(KushaShi$m12_2)
KushaShimeanmig21=mean(KushaShi$m21_2)
KushaShiNe1=KushaShi$nu1_2
KushaShiNe2=KushaShi$nu2_2

KushoShi=read.table(file="boot.Kushimoto.Shirahama.params")
mean(KushoShi$m)
median(KushoShi$m)
KushoShimeanmig=mean(KushoShi$m)
KushoShiNe1=KushoShi$nu1_1
KushoShiNe2=KushoShi$nu2_1

##Migration
allmig <- data.frame(qpcR:::cbind.na(AmaKusha$m12,AmaKusha$m21,KocKusha$m,KocShi$m12,KocShi$m21,KocKusho$m12,
             KocKusho$m21,KushaKusho$m12_2,KushaKusho$m21_2,KushaShi$m12_2,KushaShi$m21_2,KushoShi$m))
colnames(allmig)=c("AmaKusha_m12","AmaKusha_m21","KocKusha_m","KocShi_m12","KocShi_m21","KocKusho_m12","KocKusho_m21",
                   "KushaKusho_m12","KushaKusho_m21","KushaShi_m12","KushaShi_m21", "KushoShi_m")
library(reshape2)
dtam=melt(data=allmig, measure.vars=colnames(allmig))
library(dplyr)
AmaKushamig=filter(dtam, variable=="AmaKusha_m12"| variable=="AmaKusha_m21")
p1=ggplot(AmaKushamig, aes(x=variable, y=log(value)))+
  geom_boxplot()+
 # scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+
  xlab(NULL)+theme(axis.text.x =element_blank())
p1
KushaKushomig=filter(dtam, variable=="KushaKusho_m12"| variable=="KushaKusho_m21")
p2=ggplot(KushaKushomig, aes(x=factor(variable,level=c('KushaKusho_m21','KushaKusho_m12')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="Ne", labels =c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p2
KocKushomig=filter(dtam, variable=="KocKusho_m12"| variable=="KocKusho_m21")
p3=ggplot(KocKushomig, aes(x=factor(variable,level=c('KocKusho_m21','KocKusho_m12')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="Ne", labels =c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p3
KushaShimig=filter(dtam, variable=="KushaShi_m12"| variable=="KushaShi_m21")

p4=ggplot(KushaShimig, aes(x=factor(variable,level=c('KushaShi_m21','KushaShi_m12')), y=log(value)))+
  geom_boxplot()+
 # scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p4
KocShimig=filter(dtam, variable=="KocShi_m12"| variable=="KocShi_m21")
p5=ggplot(KocShimig, aes(x=factor(variable,level=c('KocShi_m21','KocShi_m12')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p5
KushoShimig=filter(dtam, variable=="KushoShi_m")
p6=ggplot(KushoShimig, aes(x=factor(variable,level=c('KushoShi_m')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p6
KocKushamig=filter(dtam, variable=="KocKusha_m")
p7=ggplot(KocKushamig, aes(x=factor(variable,level=c('KocKusha_m')), y=log(value)))+
  geom_boxplot()+
  #scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab(NULL)+xlab(NULL)+theme(axis.text.x =element_blank())
p7
library(cowplot)
plot_grid(
  p1, NULL,p7, p2, p3,NULL,p4,p5,p6,
  ncol = 3,
  label_size = 12,
  align = "v"
)

#Ne
#install.packages("qpcR")
library(qpcR)
#
Allne=cbind(AmaKushaNe1,AmaKushaNe2,KocKushaNe1,KocKushaNe2,KocShiNe1,KocShiNe2,KocKushoNe1,
            KocKushoNe2,KushaKushoNe1,KushaKushoNe2)
dta <- data.frame(qpcR:::cbind.na(AmaKushaNe1,AmaKushaNe2,KocKushaNe1,KocKushaNe2,KocShiNe1,KocShiNe2,KocKushoNe1,
                       KocKushoNe2,KushaKushoNe1,KushaKushoNe2,KushaShiNe1,KushaShiNe2,KushoShiNe1,KushoShiNe2 ))


dtane=melt(data=dta, measure.vars=colnames(dta))

ggplot(dtane, aes(x=variable, y=log(value)))+
  geom_boxplot()+
  scale_y_continuous(name="log(Ne)", labels = scales::comma)+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

Allmig=cbind(AmaKushameanmig12,AmaKushameanmig21,KocKushameanmig,KocShimeanmig12,KocShimeanmig21,KocKushomeanmig12,KocKushomeanmig21,
             KushaKushomeanmig12,KushaKushomeanmig21,KushaShimeanmig12,KushaShimeanmig21,KushoShimeanmig)
library(dplyr)
AmaKushaNes=filter(dtane, variable=="AmaKushaNe1"| variable=="AmaKushaNe2")
p1=ggplot(AmaKushaNes, aes(x=variable, y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab("Ne")+
  xlab(NULL)+theme(axis.text.x =element_blank())
p1
KushaKushoNes=filter(dtane, variable=="KushaKushoNe1"| variable=="KushaKushoNe2")
p2=ggplot(KushaKushoNes, aes(x=factor(variable,level=c('KushaKushoNe2','KushaKushoNe1')), y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Ne", labels =c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)+theme(axis.text.x =element_blank())
KocKushoNes=filter(dtane, variable=="KocKushoNe1"| variable=="KocKushoNe2")
p3=ggplot(KocKushoNes, aes(x=factor(variable,level=c('KocKushoNe2','KocKushoNe1')), y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Ne", labels =c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)+theme(axis.text.x =element_blank())
KushaShiNes=filter(dtane, variable=="KushaShiNe1"| variable=="KushaShiNe2")
p4=ggplot(KushaShiNes, aes(x=factor(variable,level=c('KushaShiNe2','KushaShiNe1')), y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)+theme(axis.text.x =element_blank())
KocShiNes=filter(dtane, variable=="KocShiNe1"| variable=="KocShiNe2")
p5=ggplot(KocShiNes, aes(x=factor(variable,level=c('KocShiNe2','KocShiNe1')), y=value))+
  geom_boxplot()+
  scale_y_continuous(name="Ne", labels = c("0","10k","20k","30k"), breaks=c(0,10000, 20000, 30000), limits=c(0,35000))+
  theme_classic(base_size = 22)+ylab("Ne")+xlab(NULL)+theme(axis.text.x =element_blank())

library(cowplot)
plot_grid(
  p1, NULL, p2, p3, p4,p5,
  ncol = 2,
  label_size = 12,
  align = "v"
)
