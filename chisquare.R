big small
##   setosa       1    49
##   versicolor  29    21
##   virginica   47     3

Sek<-data.frame(0,40)
names(Sek)<-c("clones","non.clones")

Or<-data.frame(0,38)
names(Or)<-c("clones","non.clones")

Kushima<-data.frame(2,18)
names(Kushima)<-c("clones","non.clones")
[1] 0.1111111
> .11*360
[1] 39.6

Ko<-data.frame(0,19)
names(Ko)<-c("clones","non.clones")

Shir<-data.frame(6,15)
names(Shir)<-c("clones","non.clones")

> 6/15
[1] 0.4
> .4*360
[1] 144

Kushimoto<-data.frame(6,17)
names(Kushimoto)<-c("clones","non.clones")
> 6/17
[1] 0.3529412
> 0.3529*360
[1] 127.044
Amakusa<-data.frame(6,20)
names(Amakusa)<-c("clones","non.clones")
6/20
108
newdf <- rbind(Sek, Or,Kushima,Ko, Shir,Kushimoto,Amakusa)
rownames(newdf)=c("Sekisei","Oura","Kushima","Kochi","Shirahama","Kushimoto","Amakusa")

nonmarg=data.frame(2,115)
names(nonmarg)=c("clones","non.clones")

marg=data.frame(18,52)
names(marg)=c("clones","non.clones")
newdf <- rbind(marg, nonmarg)
rownames(newdf)=c("marg","nonmarg")
#newdf$type=rownames(newdf)
chisq.test(as.matrix(newdf))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  as.matrix(newdf)
# X-squared = 23.969, df = 1, p-value = 9.789e-07
nonmarg=data.frame(2,37)
names(nonmarg)=c("clones","non.clones")

marg=data.frame(18,52)
names(marg)=c("clones","non.clones")
newdf <- rbind(marg, nonmarg)
rownames(newdf)=c("marg","nonmarg")
chisq.test(as.matrix(newdf))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  as.matrix(newdf)
# X-squared = 5.7772, df = 1, p-value = 0.01624

