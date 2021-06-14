# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

# if (length(commandArgs(trailingOnly=TRUE))<1) {
# 	options(warning.length=8000)
# 	stop("
# 	
# Summarizes bootstrap replicates for the best model.
# 
# bootRes=[filename]    results of bootstraps (likely named as [pop.contrast].winboots)
#                      
# folded=FALSE    whether the analysis was using folded SFS
# 				
# topq=0.5        top quantile to summarize. For example 0.75 means that only the best-likelihood 75% 
#                 of bootstrap replicates will be used to summarize paramter values.
#                 
# path2models=\"/projectnb/davieslab/jfifer/Japan_rad/AFS-analysis-with-moments/multimodel_inference/\"   path to the cloned repository
# 
# ")
# }

# ----------- reading input parameters

# infl=grep("bootRes=",commandArgs())
# if (length(infl)==0) { stop ("specify input file (infile=filename)\nRun script without arguments to see all options\n") }
# bootRes=sub("bootRes=","", commandArgs()[infl])
# 
# topq =grep("topq=",commandArgs())
# if(length(topq)>0) { topq=as.numeric(sub("topq=","", commandArgs()[topq])) } else { topq=0.5 }
# 
# path2models =grep("path2models=",commandArgs())
# if(length(path2models)>0) { path2models=sub("path2models=","", commandArgs()[path2models]) } else { path2models="~/AFS-analysis-with-moments/multimodel_inference/" }
# 
# if(length(grep("folded=T",commandArgs()))>0) { folded=TRUE } else { folded=FALSE }

# ----------- reading data
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimaShirahama/")
setwd("/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop2popN/")


require(ggplot2)
#bootRes="AmakusaKushima.winboots.res"
#bootRes="pop4popN.winboots.res"
bootRes="KochiKushima.winboots.res"
bootRes="KochiKushimoto.winboots.res"

bootRes="KushimaShirahama.winboots.res"

  #bootRes="KochiShirahama.winboots.res"
  path2models="/projectnb/davieslab/jfifer/Japan_rad/AFS-analysis-with-moments/multimodel_inference/"
  topq=0.5
  folded=FALSE

#system(paste("grep RESULT ", bootRes," -A 4 | grep -v Launcher | grep -E \"[0-9]|\\]\" | perl -pe 's/^100.+\\.o\\d+\\S//' | perl -pe 's/\\n//' | perl -pe 's/[\\[\\]]//g' | perl -pe 's/RESULT/\\nRESULT/g' | grep RESULT >", bootRes,".res",sep=""))

#infile=paste(bootRes,".res",sep="")
infile=bootRes
npl=read.table(infile)
pdf(paste(bootRes,"_plots.pdf",sep=""),height=3, width=8)

#------ retrieving parameter names, applying them

if (folded) { 
  pa=readLines(paste(path2models,"folded_params",sep=""))
} else {
  pa=readLines(paste(path2models,"unfolded_params",sep=""))
}
wmod=as.character(npl[1,2])
npl=npl[,-c(1:2)]
params=c(strsplit(gsub("[ \t]","",pa[grep(paste0(wmod,".py"),pa)]),split="[:,]")[[1]][-1],"theta")
names(npl)=c("id","np","ll","boot","p1","p2",params)
head(npl)
#checj this is true, if not could be some letters at the end of your value for theta. 
is.numeric(npl$theta)
#------ finding best likelihood for each bootstrap rep
maxlike=c()
for (b in 1:length(levels(npl$boot))) {
	bb=levels(npl$boot)[b]
	sub=subset(npl,boot==bb)
	maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
}
#head(maxlike)

# ---- leaving only topq % of bootstraps

hist(maxlike$ll,breaks=50,main="bootstrap likes")
abline(v=quantile(maxlike$ll,(1-topq)),col="red")
#abline(v=quantile(maxlike$ll,topq+(1-topq)/2),col="red")
maxlike=maxlike[maxlike$ll>quantile(maxlike$ll,(1-topq)),]
#hist(maxlike$ll,breaks=20)

# ---- log-transforming migration rates

migrations=grep("^m",names(maxlike))
maxlike[,migrations]=as.numeric(apply(maxlike[,migrations,drop=F],2,log,base=10))


# ----- "folding" genomic islands such that island migration is lower on average

if(length(grep("i",wmod))>0) { 
	
	mxl=maxlike
	migrations=grep("^m",names(maxlike))
	migrations.i=grep("^m.*i",names(maxlike))
	migrations=migrations[!(migrations %in% migrations.i)]

	mm=apply(mxl[,migrations],1,mean)
	mi=apply(mxl[,migrations.i],1,mean)
		
	maxlike[mm<mi,migrations]=mxl[mm<mi,migrations.i]
	maxlike[mm<mi,migrations.i]=mxl[mm<mi,migrations]
	maxlike[mm<mi,"P"]=1-maxlike[mm<mi,"P"]
}


# ---- makig a long table, defining parameter types
stack(maxlike[,7:ncol(maxlike)])
ms=stack(maxlike[,7:ncol(maxlike)])
names(ms)[2]="parameter"
ms$type="misc"
ms$type[grep("^nu",ms$parameter)]="nu"
ms$type[grep("^T",ms$parameter)]="T"
ms$type[grep("^m",ms$parameter)]="log10.migr"
ms$type[ms$parameter=="P"]="prop.islands"
ms$type[ms$parameter=="theta"]="theta"
ms$type[ms$parameter=="p_misid"]="p_misid"
ms$type=factor(ms$type,levels=c("nu","T","log10.migr","prop.islands","theta","p_misid"))
is.factor(ms$type)
# ---- plotting, saving results

pp=ggplot(ms,aes(parameter,values))+geom_boxplot()+theme(axis.text.x = element_text(angle = 45,hjust=1))+facet_wrap(~type,scale="free",nrow=1)
plot(pp)
medians=as.numeric(apply(maxlike[,7:ncol(maxlike)],2,median))


q25=apply(maxlike[,7:ncol(maxlike)],2,quantile,prob=0.25)
q75=apply(maxlike[,7:ncol(maxlike)],2,quantile,prob=0.75)

mres=data.frame(cbind(medians, q25,q75))
save(mres,maxlike,file=paste(bootRes,"_bootres.RData",sep=""))
print(medians)
dev.off()
# ---- finding represenative run (most similar to medians)

idpara=t(maxlike[,7:ncol(maxlike)])
#or idpara=t(df) if you use maxlike1 lines
cors=c()
for(i in 1:ncol(idpara)){
	cors=c(cors,cor(medians,idpara[,i]))
}
bestrun=maxlike$id[which(cors==max(cors))[1]]
system(paste("cp *_",bestrun,"_*pdf ", bootRes,"_representativeModel.pdf",sep=""))
system(paste("cp *_",bestrun,".png ", bootRes,"_representativeModel.png",sep=""))
message("representative run: ",bestrun)

###Param transform###
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KochiKushimoto/")
setwd("/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop2popN/")

load(file="pop2popN.winboots.res_bootres.RData")
#only migiration is log transformed 
mutationrate=0.018
gentime=5
antilog=10^maxlike[,15:20]
Ne=round(maxlike[,23]/(4*mutationrate),0)
antilog$nu0=round(maxlike$nu0*Ne,0)
antilog$nu1_2=round(maxlike$nu1_2*Ne,0)
antilog$nu1_3=round(maxlike$nu1_3*Ne,0)
antilog$nu2_2=round(maxlike$nu2_2*Ne,0)
antilog$nu2_3=round(maxlike$nu2_3*Ne,0)
#add the mis

antilog$m=(antilog$m*(1-maxlike[,21]))+(antilog$mi*(maxlike[,21]))
antilog$m12_3=(antilog$m12_3*(1-maxlike[,21]))+(antilog$m12_3i*(maxlike[,21]))
antilog$m21_3=(antilog$m21_3*(1-maxlike[,21]))+(antilog$m21_3i*(maxlike[,21]))



antilog$m=round(antilog$m/(Ne*2),6)
antilog$m12_3=round(antilog$m12_3/(Ne*2),6)
antilog$m21_3=round(antilog$m21_3/(Ne*2),6)

antilog$T2=round(gentime*maxlike$T3*2*Ne,0)
antilog$T1=round(gentime*maxlike$T3*2*Ne,0)+round(gentime*maxlike$T2*2*Ne,0)
antilog$T0=round(gentime*maxlike$T3*2*Ne,0)+round(gentime*maxlike$T2*2*Ne,0)+round(gentime*maxlike$T1*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.pop2popN.out.params")
#######
setwd("/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop4popN/")

load(file="pop4popN.winboots.res_bootres.RData")
antilog=10^maxlike[,13:16]
Ne=round(maxlike$theta/(4*mutationrate),0)

antilog$nu1_2=round(maxlike$nu1_2*Ne,0)
antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_2=round(maxlike$nu2_2*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)

antilog$m12_2=(antilog$m12_2*(1-maxlike$P))+(antilog$m12_2i*(maxlike$P))
antilog$m21_2=(antilog$m21_2*(1-maxlike$P))+(antilog$m21_2i*(maxlike$P))


antilog$m12_2=round(antilog$m12_2/(Ne*2),6)
antilog$m21_2=round(antilog$m21_2/(Ne*2),6)

#no migration in epoch 1
antilog$T1=round(gentime*maxlike$T*2*Ne,0)
antilog$T0=round(gentime*maxlike$T*2*Ne,0)+round(gentime*maxlike$T0*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.pop4popN.out.params")


###
#######
setwd("/projectnb/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop2pop4/")

load(file="pop2pop4.winboots.res_bootres.RData")
antilog=10^maxlike[,15:20]
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$nu0=round(maxlike$nu0*Ne,0)
antilog$nu1_2=round(maxlike$nu1_2*Ne,0)
antilog$nu1_3=round(maxlike$nu1_3*Ne,0)
antilog$nu2_2=round(maxlike$nu2_2*Ne,0)
antilog$nu2_3=round(maxlike$nu2_3*Ne,0)

antilog$m=(antilog$m*(1-maxlike$P))+(antilog$mi*(maxlike$P))
antilog$m12_3=(antilog$m12_3*(1-maxlike$P))+(antilog$m12_3i*(maxlike$P))
antilog$m21_3=(antilog$m21_3*(1-maxlike$P))+(antilog$m21_3i*(maxlike$P))

antilog$m=round(antilog$m/(Ne*2),6)
antilog$m12_3=round(antilog$m12_3/(Ne*2),6)
antilog$m21_3=round(antilog$m21_3/(Ne*2),6)

antilog$T2=round(gentime*maxlike$T3*2*Ne,0)
antilog$T1=round(gentime*maxlike$T3*2*Ne,0)+round(gentime*maxlike$T2*2*Ne,0)
antilog$T0=round(gentime*maxlike$T3*2*Ne,0)+round(gentime*maxlike$T2*2*Ne,0)+round(gentime*maxlike$T1*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.pop2pop4.out.params")
####
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.AmakusaKushima/")
load(file="AmakusaKushima.winboots.res_bootres.RData")
antilog=10^maxlike[,10:11]
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)
antilog$m12=round(antilog$m12/(Ne*2),6)
antilog$m21=round(antilog$m21/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.AmakusaKushima.out.params")
####
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KochiKushima/")
load(file="KochiKushima.winboots.res_bootres.RData")
antilog=as.data.frame(10^maxlike[,11])
colnames(antilog)="m"
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)
antilog$m=round(antilog$m/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T0*2*Ne,0)
antilog$T1=round(gentime*maxlike$T0*2*Ne,0)+round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.KochiKushima.out.params")
###
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KochiKushimoto/")
load(file="KochiKushimoto.winboots.res_bootres.RData")
antilog=10^maxlike[,10:11]
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)
antilog$m12=round(antilog$m12/(Ne*2),6)
antilog$m21=round(antilog$m21/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.KochiKushimoto.out.params")
###
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KochiShirahama/")
load(file="KochiShirahama.winboots.res_bootres.RData")
antilog=10^maxlike[,10:11]
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)
antilog$m12=round(antilog$m12/(Ne*2),6)
antilog$m21=round(antilog$m21/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.KochiShirahama.out.params")
###
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimaKushimoto/")
load(file="KushimaKushimoto.winboots.res_bootres.RData")
antilog=10^maxlike[,12:13]
Ne=round(maxlike$theta/(4*mutationrate),0)
antilog$n0=round(maxlike$nu0*Ne,0)
antilog$nu1_2=round(maxlike$nu1_2*Ne,0)
antilog$nu2_2=round(maxlike$nu2_2*Ne,0)
antilog$m12_2=round(antilog$m12_2/(Ne*2),6)
antilog$m21_2=round(antilog$m21_2/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T0*2*Ne,0)
antilog$T1=round(gentime*maxlike$T0*2*Ne,0)+round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.KushimaKushimoto.out.params")

### NEED TO COME BACK TO THIS ONE 
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimaShirahama/")
load(file="KushimaShirahama.winboots.res_bootres.RData")
antilog=10^maxlike[,11:18]
Ne=round(maxlike$theta/(4*mutationrate),0)

antilog$nu1_2=round(maxlike$nu1_1*Ne,0)
antilog$nu2_2=round(maxlike$nu2_1*Ne,0)

antilog$m12_1=(antilog$m12_1*(1-maxlike$P))+(antilog$m12_1i*(maxlike$P))
antilog$m21_1=(antilog$m21_1*(1-maxlike$P))+(antilog$m21_1i*(maxlike$P))
antilog$m12_2=(antilog$m12_2*(1-maxlike$P))+(antilog$m12_2i*(maxlike$P))
antilog$m21_2=(antilog$m21_2*(1-maxlike$P))+(antilog$m21_2i*(maxlike$P))

antilog$m12_2=round(antilog$m12_2/(Ne*2),6)
antilog$m21_2=round(antilog$m21_2/(Ne*2),6)

antilog$T0=round(gentime*maxlike$T*2*Ne,0)
antilog$T1=round(gentime*maxlike$T0*2*Ne,0)+round(gentime*maxlike$T*2*Ne,0)
antilog$Ne=Ne
antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.KushimaShirahama.out.params")
######
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimotoShirahama/")
load(file="KushimotoShirahama.winboots.res_bootres.RData")
antilog=as.data.frame(10^maxlike[,11])
colnames(antilog)="m"
Ne=round(maxlike$theta/(4*mutationrate),0)

antilog$nu1_1=round(maxlike$nu1_1*Ne,0)
antilog$nu2_1=round(maxlike$nu2_1*Ne,0)

antilog$m=round(antilog$m/(Ne*2),6)


antilog$T0=round(gentime*maxlike$T*2*Ne,0)
antilog$T1=round(gentime*maxlike$T0*2*Ne,0)+round(gentime*maxlike$T*2*Ne,0)

antilog$ll=maxlike$ll
antilog$p1=maxlike$p1
antilog$p2=maxlike$p2
write.table(antilog,file="boot.Kushimoto.Shirahama.params")
