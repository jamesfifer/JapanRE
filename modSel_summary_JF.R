# read about multimodel inference here: 
# https://pdfs.semanticscholar.org/a696/9a3b5720162eaa75deec3a607a9746dae95e.pdf

# if (length(commandArgs(trailingOnly=TRUE))<1) {
# 	options(warning.length=8000)
# 	stop("
# 	
# Summarizes results of the model selection run, 
# then writes a list of commands for bootstrapping the winning model to the file [pop.contrast].winboots.runs
# 
# Model selection is assumed to have been run on 10+ boostrapped replicates, each model with 3+ random starts.
# 
# Arguments:
# 
# modselResult=[filename]          file containing results of the model selection run. 
# 
# The following parameters are about writing commands for bootstrapping the winning model
# (NOT about the completed run for model selection):
# 
# nboots=100                       number of boostrap replicates to run for the winning model
# nreps=6                          number of random restarts for each boostrap replicate
# folded=FALSE                     whether the analysis is using folded SFS
# args=[list of arguments]           names of pop1 and pop2, projection for pop1, projection for pop2, mutation rate (per genotyped portion of the 
#                                    genome per generation), generation time in thousands of years. Population names can be anything. 
#                                    For ANGSD-derived SFS, projections should be 0.8*2N for each population (ronded to integer); 
#                                    in the case corresponding to the default setting each population was represented by 10 individuals.
#                                    Example: args=\"p1 p2 16 16 0.02 0.005\"
# 
# ")
# }
# 
# infl=grep("modselResult=",commandArgs())
# if (length(infl)==0) { stop ("specify module selection results file (modselResult=filename)\nRun script without arguments to see all options\n") }
# modselResult=sub("modselResult=","", commandArgs()[infl])
# 
# args=grep("args=",commandArgs())
# if (length(args)==0) { stop ("specify arguments for SFS model runs\nRun script without options to see details\n") }
# args =sub("args=","", commandArgs()[args])
# 
# if(length(grep("folded=T",commandArgs()))>0) { folded=TRUE } else { folded=FALSE }
# 
# nreps =grep("nreps=",commandArgs())
# if(length(nreps)>0) { nreps=as.numeric(sub("nreps=","", commandArgs()[nreps])) } else { nreps=6 }
# 
# nboots =grep("nboots=",commandArgs())
# if(length(nboots)>0) { nreps=as.numeric(sub("nboots =","", commandArgs()[nboots])) } else { nboots=100 }
# 
# if(length(grep("folded=T",commandArgs()))>0) { folded=TRUE } else { folded=FALSE }
#comargs = commandArgs(trailingOnly=TRUE)
setwd("/projectnb2/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop2pop4/")
setwd("/projectnb2/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop2popN/")
setwd("/projectnb2/davieslab/jfifer/Japan_rad/admixed_pops_analysis/4pop_analysis/Demographic_Analysis_II/Unlinked_analyses/boot.pop4popN/")

setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.AmakusaKushima")
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimotoShirahama")
setwd("/projectnb/davieslab/jfifer/Japan_rad/Demographic_analysis_II/North_Only/Unlinked_Analyses/boot.KushimaShirahama")

modselResult="KushimaShirahama.res"
#modselResult=comargs[1] 
path2models="/projectnb/davieslab/jfifer/Japan_rad/AFS-analysis-with-moments/multimodel_inference/"
 folded=FALSE
 #args="pop2 popN 33 150 0.018 0.005"
 #args="pop2 pop4 33 68 0.018 0.005"
#args="pop4 popN 68 150 0.018 0.005"
 #args="Amakusa Kushima 29 32 0.018 0.005"
 #args="Kochi Shirahama 30 29 0.018 0.005"
 args= "Kochi Kushima 30 32 0.018 0.005"
 #args="Kochi Kushimoto 30 27 0.018 0.005"
 #args="Kushima Kushimoto 32 27 0.018 0.005"
 args="Kushima Shirahama 32 29 0.018 0.005"
# args="Kushimoto Shirahama 27 29 0.018 0.005"
 #args=comargs[2]
 nreps=6
 nboots=100

#system(paste("grep RESULT ", modselResult," -A 4 | grep -v Launcher | grep -E \"[0-9]|\\]\" | perl -pe 's/^100.+\\.o\\d+\\S//' | perl -pe 's/\\n//' | perl -pe 's/[\\[\\]]//g' | perl -pe 's/RESULT/\\nRESULT/g' | grep RESULT >", modselResult,".res",sep=""))
#infile=paste(modselResult,".res",sep="")
infile=modselResult

library(ggplot2)

system(paste("cut -f 2,3,4,5,6 -d ' ' ",infile," > ",infile,".likes",sep=""))
npl=read.table(paste(infile,".likes",sep=""))
#read.table(file="KochiKushima.res.likes")
system(paste("rm ",infile,".likes",sep=""))

names(npl)=c("model","id","npara","ll","boot")
contrast=sub("_.+","", npl$boot[1])

#head(npl)
#npl=npl[grep(infile,npl$boot),]
aics=list()
for (b in 1:length(unique(npl$boot))) {
	bb=levels(npl$boot)[b]
	nplb=subset(npl,boot==bb)
	maxlike=c();nmod=c()
	for (m in unique(nplb$model)) {
		sub=subset(nplb,model==m)
		nmod=c(nmod,nrow(sub))
		maxlike=data.frame(rbind(maxlike,sub[sub$ll==max(sub$ll),]))
	}
	npara=maxlike$npara
	likes=maxlike$ll
	aic=2*npara-2*likes
	aicc=data.frame(cbind(model=as.character(unique(nplb$model))))
	aicc$aic=aic
	aicc$nmod=nmod
	aicc$boot=bb
	aics[[b]]=aicc
}
awt=data.frame(do.call(rbind,aics))
models=unique(awt$model)
med=c()
for (m in models) {
	ss=subset(awt,model==m)
	med=c(med,median(ss$aic))
}
modmed=data.frame(cbind(mod=as.character(models)))
modmed$med=med
modmed=modmed[order(med),]
modmed$mod=factor(modmed$mod,levels=modmed$mod)
awt$model=factor(awt$model,levels=modmed$mod)

pdf(width=2.2,height=2,file=paste(contrast,"_modsel_top10medians.pdf",sep=""))
pp=ggplot(modmed[1:10,],aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
#ggplot(modmed,aes(mod, med))+geom_point()+theme(axis.text.x = element_text(angle = 45,hjust=1))
plot(pp)
dev.off()
pdf(width=15,height=5,file=paste(contrast,"_modsel_allBoxplots.pdf",sep=""))
pp=ggplot(awt,aes(model,aic))+geom_boxplot()+scale_y_log10()+theme(axis.text.x = element_text(angle = 45,hjust=1))
plot(pp)
dev.off()

# ----- extracting name and parameters of the winning model, writing them to a file
#For Kushima Shirahama, two are basically the same (sc1i, sc11sm) so choosing sc11sm
#winner=as.character(modmed[2,1])

winner=as.character(modmed[1,1])
winner
npl0=subset(npl,model==winner)
npl0=npl0[which(npl0$ll==max(npl0$ll)),]
npl0
system(paste("grep -a \"",npl0$id," \" ", infile," > ",infile,".winmod",sep=""))

system(paste("rm ",infile,sep=""))

npl0=read.table(paste(infile,".winmod",sep=""))
params=as.vector(npl0[1,c(9:(ncol(npl0)-1))])
write.table(params,file=paste(contrast,".",winner,sep=""),quote=F,col.names=F,row.names=F)

# print(as.character(modmed[1,1]),quote=F)

# ------ writing commands to bootstrap the winning model

# if (folded) { 
# 	mods=paste("fold_",winner,sep="")
# 	} else {
# 	mods=winner 		
# 	}
# 
# args2=c()
# for (b in 1:nboots) {
# 	for (n in 1:nreps) {
# 		for (m in mods) {
# 			bname=paste(contrast,"_",b,".sfs",sep="")
# 			args2=c(args2,paste("sleep ",n," && ",m,".py ",bname," ",args," ",paste(contrast,".",winner,sep="")," >>",contrast,".winboots",sep=""))
# 		}
# 	}
# }
# 
# writeLines(args2,paste(contrast,".winboots.runs",sep=""))



