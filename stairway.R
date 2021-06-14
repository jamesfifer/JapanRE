setwd("~/BOSTON/Davies/Range expansion/Code/stairwayplot")
popN=read.table(file="popN.final.summary",header=T)
popN$population="popN"
pop2=read.table(file="two-epoch.final.summary",header=T)
pop2$population="pop2"
pop4=read.table(file="pop4.final.summary",header=T)
pop4$population="pop4"

all=rbind(popN,pop2,pop4)
library(ggplot2)
ggplot(data=all, aes(x=year, y=Ne_median, color=population)) +geom_line(size=1.5)+
geom_ribbon(aes(ymin=all$Ne_12.5.,ymax=all$Ne_87.5., fill = factor(population)), alpha = 0.3,colour = NA,show_guide = FALSE)+
 theme_bw(base_size = 22)+ scale_y_continuous(name="Ne", labels = scales::comma)+#scale_x_continuous(name="Years", labels=c("0","100k","200k",
  scale_x_continuous(name="Years",limits =  c(0,6e+05))+                                                                                                                         # "300k"),limits =  c(0,3e+05))+ 
  scale_fill_manual(values=c("#caff70", "#006400", "#ff0088"))+
  scale_color_manual(values=c("#caff70", "#006400", "#ff0088"))
  






