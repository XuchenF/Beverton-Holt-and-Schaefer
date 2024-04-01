k<-100000
M=1.19
h<-0.36
hh<-0.99

b<-(h-0.2)/(0.2*k-0.2*h*k)
a<-(k*(1-exp(-1*M))*(1+b*k))/k

bb<-(hh-0.2)/(0.2*k-0.2*hh*k)
aa<-(k*(1-exp(-1*M))*(1+bb*k))/k

#import a col from the excel spreedsheet into the R
B<-read.csv("Biomass-simulation.csv")
#Using Beverton-Holt driven production model to calculate corresponding ideal catch
##All the parameters are based on the Yellow Sea anchovy information
#h=0.36
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

#hh=0.99
C<-data.frame(matrix(ncol=6, nrow=50))
F<-data.frame(matrix(ncol=6, nrow=50))
for(j in 1:6){
for(i in 1:50){
if (B[i+1,j] > aa*B[i,j]/(1+bb*B[i,j])){
F[i,j]=(-1)*log((B[i+1,j]-(aa*B[i,j]/(1+bb*B[i,j])))/B[i,j])-M
C[i,j]=(F[i,j]/(F[i,j]+M))*B[i,j]*(1-exp(-(F[i,j]+M)))
}
else
{
C[i,j]=((log(10*B[i,j])-M)/log(10*B[i,j]))*B[i,j]+(aa*B[i,j]/(1+bb*B[i,j]))-B[i+1,j]
}
}
}


