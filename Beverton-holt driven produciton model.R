k<-100000
M=0.9
h<-0.66

b<-(h-0.2)/(0.2*k-0.2*h*k)
a<-(k*(1-exp(-1*M))*(1+b*k))/k


#import a col from the excel spreedsheet into the R
B<-read.csv("Biomass-simulation.csv")
#Using Beverton-Holt driven production model to calculate corresponding ideal catch
##All the parameters are based on the Yellow Sea anchovy information
C<-data.frame(matrix(ncol=6, nrow=50))
F<-data.frame(matrix(ncol=6, nrow=50))
for(j in 1:6){
for(i in 1:50){
if (B[i+1,j] > a*B[i,j]/(1+b*B[i,j])){
F[i,j]=(-1)*log((B[i+1,j]-(a*B[i,j]/(1+b*B[i,j])))/B[i,j])-M
C[i,j]=(F[i,j]/(F[i,j]+M))*B[i,j]*(1-exp(-(F[i,j]+M)))
}
else
{
C[i,j]=((log(10*B[i,j])-M)/log(10*B[i,j]))*B[i,j]+(a*B[i,j]/(1+b*B[i,j]))-B[i+1,j]
}
}
}


Binit<-c((max(C[,2])/1.76+4*max(C[,2])/0.98)/2*0.6, Binit<-(max(C[,2])/1.76+4*max(C[,2])/0.98)/2*0.8)
for(i in 1:length(C[,1])){
if(i-1<1){
y=function(x){x/(x+0.9)*Binit*(1-exp(-1*(x+0.9)))-C[i,2]}
so<-uniroot(y, interval = c(0, 17))
F[i,1]<-so$root
B[i,1]<-a*Binit/(1+b*Binit)+Binit*exp(-1*(F[i,1]+0.9))
}
else
{
y=function(x){x/(x+0.9)*B[i-1,1]*(1-exp(-1*(x+0.9)))-C[i,2]}
so<-uniroot(y, interval = c(0, 17))
F[i,1]<-so$root
B[i,1]<-a*Binit/(1+b*B[i-1,1])+B[i-1,1]*exp(-1*(F[i,1]+0.9))
}
}
