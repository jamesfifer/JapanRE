#!/bin/env Rscript

#BiocManager::install("snpStats")
library(snpStats)
#devtools::install_github("BenjaminPeter/rangeExpansion", ref="package")
library(rangeExpansion) 
library(sp)
region <- list(NULL)
ploidy <- 2

#Acropora.snp="sfsSites.snapp"
#Acropora.coord="ACTUALFINALCOORD.csv" #without kushimoto (pop4)
#Acropora.snp="sfsSites_nokush.snapp"
#Acropora.coord="ACTUALFINALCOORD_5reg_nokush.csv"
#without paralogs, without fixed ancestrals
Acropora.snp="fixed.anc.rm.snapp"
Acropora.coord="ACTUALFINALCOORD.csv"

#next time try filtering one snp per rad tag


raw.data <- load.data.snapp(Acropora.snp,Acropora.coord,sep=',', ploidy=ploidy)
#

print(raw.data)

pop <- make.pop(raw.data, ploidy) 
print(pop)
psi <- get.all.psi(pop)
res <- run.regions(region=region, pop=pop, psi=psi)
#We calculated the directionality index, ψ, for all population pairs using the get.all.psi function. To determine significance, we calculated the 
#standard error of the upper triangle of the pairwise ψ matrix excluding the diagonal, thereby allowing us to calculate the Z-score for each population. 
#For each region of interest, we plotted data for each pair of populations where the absolute Z-score was greater than 5 and visually assessed the 
#geographic patterns of source and sink populations.
save.image(file='slatkin3.RData')
summary(res)
#get pairwise estimates of psi
psi

#export plot
pdf(file="NAME.all.pdf")

plot(res)
dev.off() 



