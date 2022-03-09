data1<- read.csv("BTC-USD.csv")
data2<- read.csv("ETH-USD.csv")


head(data1)
head(data2)

View(data1)
View(data2) 

colnames(data1)[1] <- "Date"
colnames(data1)[7] <- "Change"
colnames(data2)[1] <- "Date"
colnames(data2)[7] <- "Change"

data1$Price<- sub(",","", data1$Price)
data1$Price = as.double(ingl.price)
#spcjt.price <- data2$Price
#summary(ingl.price)
#summary(spcjt.price)
summary(data2)
summary(ingl.price)

mean(ingl.price)

ts.plot(ingl.price)

#############data1 plots#################################
data1$Date = as.Date(data1$Date, format="%Y-%m-%d")
data1<-data1[order(as.Date(data1$Date, format="%Y-%m-%d")),]
plot(data1$Date, data1$Close, type = 'l')

log1 = numeric(length = dim(data1)[1] -1)
for(i in 2:dim(data1)[1]){
  log1[i-1]  = log(data1$Close[i]) - log(data1$Close[i-1])
  
}
plot(log1, type = 'l')


##############data2 plots#################

data2$Date = as.Date(data2$Date, format="%Y-%m-%d")
data2<-data2[order(as.Date(data2$Date, format="%Y-%m-%d")),]
plot(data2$Date, data2$Close, type = 'l')
# log return
log2 = numeric(length = dim(data2)[1] -1)
for(i in 2:dim(data2)[1]){
  log2[i-1]  = log(data2$Close[i]) - log(data2$Close[i-1])
  
}
plot(log2, type = 'l')

###########################################
install.packages("cubfits")
library(cubfits)  # for asymmetric distribution
#plot(seq(-2,2 ,by = 0.1 ),dasla(seq(-2,2 ,by = 0.1 ), kappa = .5 , sigma =1 ), type = 'l')
library(mvtnorm)
rasymlp <- function(n=1, del.l, del.u){
    return (rasla(n,kappa = sqrt(del.l/del.u) , sigma = sqrt(1/(del.l*del.u))))
}
dasymlp <- function( x,del.l,del.u){
   if(x <=0){
      return ((del.l/(del.l + del.u))*exp(x/del.l))
   }
  else {
    return (1 - ((del.u/(del.l + del.u))*exp(-1*x/del.u)) )
  }
}

## Doubt - giving streamline structure when lower del.l and del.u ##

x1 = rasymlp(1,1-del.l,1- del.u) + r
x2 = rasymlp(1, 1- del.l, 1-del.u) +r
r = asymlp(1,del.l, del.u)
plot(x1,x2) 
##
cop <- function(w1, w2, del.l, del.u, rho){
q = c(qnorm(dasymlp(w1 , 1- del.l, 1-del.u)) , qnorm(dasymlp(w2 , 1- del.l, 1-del.u)))
ret = pmvnorm(upper = q, mean = c(0,0), sigma  = matrix(c(1,rho,rho,1), nrow =2) )
 return (ret[1])
} 

# estimating c(t,t) using monte carlo method
del.l =  0.7
del.u  = 0.2
#define x1 x2
#x1 = data
#x2 = data

t  = seq(0.1,0.99, by = 0.01)
x.u = numeric(length(t))
x.l = numeric(length(t))
for(i in 1:length(t)){
jointp = numeric(length = 1e2)
for( j in 1:1e2){
  r = rasymlp(1,del.l, del.u)
  jointp[j] =   cop(t[i] -r , t[i]-r , del.l , del.u, rho= 0.5) #C(t,t)
}
C.t.t = mean(jointp)
x.l[i] = C.t.t/t[i]
x.u[i] = (1- 2*t[i] + C.t.t)/(1-t[i])
}

 plot(t, x.l, type ='l')
plot(t,x.u, type = 'l')
