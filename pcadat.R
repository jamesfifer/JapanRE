setwd("/projectnb/davieslab/jfifer/Japan_rad/Bayescan")
#install.packages("pcadapt")

library(pcadapt)
#converted vcfs to bed using plink --vcf myresult13.vcf --make-bed --allow-extra-chr --out myresult13
#filtered is without the weird north sites

#infile="your_filtered_snpsnorth4.bed"
#infile="your_filtered_snpsnorth2.bed"
#infile="your_filtered_snps13.bed"
infile="myresult24.bed"


filename <- read.pcadapt(input = infile, type = "bed")
x <- pcadapt(input = filename, K = 20)

plot(x, option = "screeplot")
#choose left of straight line,  k=4 here
#can also based on admixture which would be 2

#check bspops files

#north4
poplist.names <- c(rep("North", 92),rep("Pop4", 43))

#North2
# With names
poplist.names <- c(rep("North", 92),rep("Pop2", 21))
print(poplist.names)

#13
poplist.names <- c(rep("Pop1", 27),rep("Pop3", 53))
#24
poplist.names <- c(rep("Pop2", 21),rep("Pop4", 43))

plot(x, option = "scores", pop = poplist.names)
#check this to see if K should be greater than 2, if there are still separate clusters should increase K
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names,plt.pkg = "plotly")
#K=2 for all comparisons
#if outliers
#vcftools --remove-indv ind11 --remove-indv ind13 --vcf myresultnorth2.vcf --recode --out your_filtered_snpsnorth2.vcf
# plink --vcf your_filtered_snpsnorth2.vcf.recode.vcf --make-bed --allow-extra-chr --out your_filtered_snpsnorth2

x <- pcadapt(filename, K = 2)
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

#lots of other transformations you can do with p values, I did q values for bayescan with FDR of 0.05
#so i did the same here to keep it consistent 
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
outliers_qvals=qval[outliers]
length(outliers)
outliers

#bim=read.table(file="your_filtered_snpsnorth4.bim")
#bim=read.table(file="your_filtered_snpsnorth2.bim")
#bim=read.table(file="your_filtered_snps13.bim")
bim=read.table(file="myresult24.bim")

outliers.bim <- bim[outliers,]
positions.bim=outliers.bim[c(1,4)]
positions.bim$V1<- sub('(?=^[0-9])', 'chr', positions.bim$V1, perl = T)
positions.bim$V5=qval[outliers]

#write.table(positions.bim, file="pcadapt_outliers_north2.txt",col.names=FALSE)
#write.table(positions.bim, file="pcadapt_outliers_north4.txt",col.names=FALSE)
#write.table(positions.bim, file="pcadapt_outliers_13.txt",col.names=FALSE)
write.table(positions.bim, file="pcadapt_outliers_24.txt",col.names=FALSE)
