library(reshape2)
library(ggplot2)
library(ggpubr)

zerop<-matrix(, nrow=50, ncol=100)
for(i in 1:50){
  zerop[i,]=C[i]*exp(rnorm(100, mean=0, sd=0.05))
} #Before execute this step, please sure we have calculated the catch by applying production funtion to the simulated biomass
zeropp<-melt(zerop, id=1)
zeropp<-as.data.frame(zeropp)
zeropp<-zeropp[, c("Var2","Var1","value")]
colnames(zeropp)<-c("Stock","yr","ct")
#----------------------------------------------------------------------
onep<-matrix(, nrow=50, ncol=100)
for(i in 1:50){
  onep[i,]=C[i]*exp(rnorm(100, mean=0, sd=0.15))
}
onepp<-melt(onep, id=1)
onepp<-as.data.frame(onepp)
onepp<-onepp[, c("Var2","Var1","value")]
colnames(onepp)<-c("Stock","yr","ct")
#----------------------------------------------------------------------
secondp<-matrix(, nrow=50, ncol=100)
for(i in 1:50){
  secondp[i,]=C[i]*exp(rnorm(100, mean=0, sd=0.25))
}
secondpp<-melt(secondp, id=1)
secondpp<-as.data.frame(secondpp)
secondpp<-secondpp[, c("Var2","Var1","value")]
colnames(secondpp)<-c("Stock","yr","ct")
#----------------------------------------------------------------------
thirdp<-matrix(, nrow=50, ncol=100)
for(i in 1:50){
  thirdp[i,]=C[i]*exp(rnorm(100, mean=0, sd=0.35))
}
thirdpp<-melt(thirdp, id=1)
thirdpp<-as.data.frame(thirdpp)
thirdpp<-thirdpp[, c("Var2","Var1","value")]
colnames(thirdpp)<-c("Stock","yr","ct")
#----------------------------------------------------------------------
forthp<-matrix(, nrow=50, ncol=100)
for(i in 1:50){
  forthp[i,]=C[i]*exp(rnorm(100, mean=0, sd=0.45))
}
forthpp<-melt(forthp, id=1)
forthpp<-as.data.frame(forthpp)
forthpp<-forthpp[, c("Var2","Var1","value")]
colnames(forthpp)<-c("Stock","yr","ct")



yr<-c(1:length(C))
data<-data.frame(yr, C)

p1<-ggplot()+geom_line(zeropp, mapping=aes(x=yr, y=ct, group=Stock), color="#E0E0E0")+geom_line(data, mapping=aes(x=yr, y=C),color="green")+theme_bw()+scale_x_continuous(name = "year", limits = c(1, 50),breaks=c(1,10,20,30,40,50))+ylim(0, 110000)+theme(legend.position = 'none')+ylab("catch")
p2<-ggplot()+geom_line(onepp, mapping=aes(x=yr, y=ct, group=Stock), color="#E0E0E0")+geom_line(data, mapping=aes(x=yr, y=C),color="green")+theme_bw()+scale_x_continuous(name = "year", limits = c(1, 50),breaks=c(1,10,20,30,40,50))+ylim(0, 110000)+theme(legend.position = 'none')+ylab("catch")
p3<-ggplot()+geom_line(secondpp, mapping=aes(x=yr, y=ct, group=Stock), color="#E0E0E0")+geom_line(data, mapping=aes(x=yr, y=C),color="green")+theme_bw()+scale_x_continuous(name = "year", limits = c(1, 50),breaks=c(1,10,20,30,40,50))+ylim(0, 110000)+theme(legend.position = 'none')+ylab("catch")
p4<-ggplot()+geom_line(thirdpp, mapping=aes(x=yr, y=ct, group=Stock), color="#E0E0E0")+geom_line(data, mapping=aes(x=yr, y=C),color="green")+theme_bw()+scale_x_continuous(name = "year", limits = c(1, 50),breaks=c(1,10,20,30,40,50))+ylim(0, 110000)+theme(legend.position = 'none')+ylab("catch")
p5<-ggplot()+geom_line(forthpp, mapping=aes(x=yr, y=ct, group=Stock), color="#E0E0E0")+geom_line(data, mapping=aes(x=yr, y=C),color="green")+theme_bw()+scale_x_continuous(name = "year", limits = c(1, 50),breaks=c(1,10,20,30,40,50))+ylim(0, 110000)+theme(legend.position = 'none')+ylab("catch")
ggarrange(p1, p2, p3, p4, p5, labels=NULL, nrow=2, ncol=3)
