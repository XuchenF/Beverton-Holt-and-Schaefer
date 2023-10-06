#import a col from the excel spreedsheet into the R
B<-read.csv("Biomass-simulation.csv")
B<-B[,1]
#Using Beverton-Holt driven production model to calculate corresponding ideal catch
##All the parameters are based on the Yellow Sea anchovy information
C<-vector( "numeric" , length(B)-1 )
for(i in 1:length(C)){
if (B[i+1] > 1.566*B[i]/(1+0.0000125*B[i])){
F[i]=(-1)*log((B[i+1]-(1.566*B[i]/(1+0.0000125*B[i])))/B[i])-1.19
C[i]=(F[i]/(F[i]+1.19))*B[i]*(1-exp(-(F[i]+1.19)))
}
else
{
C[i]=((log(10*B[i])-1.19)/log(10*B[i]))*B[i]+(1.566*B[i]/(1+0.0000125*B[i]))-B[i+1]
}
}
