data1<- read.csv("btc.csv")
data2<- read.csv("eth.csv")

head(data1)
head(data2)
View(data1)
View(data2)

data1$Price<- sub(",","", data1$Price)
btc.price = as.double(data1$Price)
#scale(btc.price, center = TRUE, scale = TRUE)
data2$Price<- sub(",","", data2$Price)
eth.price = as.double(data2$Price)
#scale(eth.price, center = TRUE, scale = TRUE)
plot(scale(btc.price, center = TRUE),scale(eth.price, center = TRUE, scale = TRUE))

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

##################libraries############
install.packages("cubfits")
install.packages("mvtnorm")
library(cubfits)  # for asymmetric distribution
#plot(seq(-2,2 ,by = 0.1 ),dasla(seq(-2,2 ,by = 0.1 ), kappa = .5 , sigma =1 ), type = 'l')
library(mvtnorm)

###########plot check##########################
x = seq(-4,4,by = 0.1)
plot(x , dasla(x, kappa = 1 , sigma = sqrt(2)), type = 'l')   # here sigma = sqrt(2*del.l*del.u) by definition of function in library

############Asymmetreic#Laplace###########################
#sampling
rasymlp <- function(n=1, del.l, del.u){
  sample = rasla(n, theta = 0,kappa = sqrt(del.l/del.u) , sigma = sqrt(2*del.l*del.u))
  return (sample)
}
#distribution
# dasymlp <- function( x,del.l,del.u){

#   solve <- function(x, del.l,del.u){
#   if(x <=0){
#     dist = (exp(x/del.l))*(del.l/(del.l + del.u))
#     return (dist)
#   }
#   else {
#     dist = (exp(-x/del.u))*(del.u/(del.l+del.u))
#     dist = 1 - dist
#     return (dist)
#   }
#   }
#   return(sapply(x, solve, del.l = del.l ,del.u = del.u))
# } 
dasymlp <- function(x,del.l,del.u){
  coeff <- del.l*as.numeric(x<=0) + del.u*as.numeric( x>0)
  
  out <- (coeff*exp(-abs(x)/coeff))/(del.l+del.u)
  out <- ifelse(x<=0 , out, 1-out)
  return(out)
}
#print(dasymlp(seq(-2,2,0.2) , 0.4,0.5))
#density
denasymlp <- function(x,del.l, del.u){
  coeff <- del.l*as.numeric(x<=0) + del.u*as.numeric(x>0)
  out <-  exp(-abs(x)/coeff)/(del.l + del.u)
  return (out)
  
}
#denasymlp(seq(-2,2,0.2), 0.4, 0.5)
###########plot check###################
#set.seed(121)
#del.l = 0.7
#del.u = 0.8
#r = rasymlp(1000,del.l, del.u)
#x1 = rasymlp(1000,1-del.l,1- del.u) + r
#x2 = rasymlp(1000, 1- del.l, 1-del.u) +r

#plot(x1,x2) 
############################
k1 = function(del.l , del.u){
  k = (del.l^3)/((del.l + del.u)*(2*del.l -1 )*(1 + del.l - del.u))
  return (k)
}
#  k1(4,3)

k2  = function(del.l,del.u){
  k = ((del.l - 1)^3)/((2*del.l -1)*(del.l -del.u -1)*(2 - del.l -del.u))
  return (k)
}
# k2(4,3)

k3 = function(del.l, del.u){
  k= (del.u^3)/((del.l + del.u)*(2*del.u -1)*(del.l - del.u -1))
  return (k)
}
#k3(2,3)
k4 = function(del.l, del.u){
  k = ((del.u -1)^3)/((2*del.u -1)*(1 + del.l - del.u)*(2- del.l - del.u))  # last breacket might be (2-del.u-del.u) as per paper
  return (k)
}

k5 = function(del.l, del.u){
  out  =  2/((1+2*del.u)*(2*del.u -3))
  return (out)
}
k6 = function(del.l, del.u){
  out  =   ((-12*del.u + 12*del.u*del.u  -5))/((1+2*del.u)*(1+2*del.u)*(2*del.u - 3)*(2*del.u - 3))
  return (out)
}
k7 = function(del.l, del.u){
  out  =  ((4*(del.u**3))/((1+2*del.u)*(1+2*del.u)*(2*del.u -1)))
  return (out)
}
k8 = function(del.l, del.u){
  out  =  (4*((del.u -1)**3))/((2*del.u-3)*(2*del.u-3)*(2*del.u -1))
  return (out)
}
k9 = function(del.l, del.u){
  out  = ((4*(del.l**3))/((1+2*del.l)*(1+2*del.l)*(2*del.l -1)))
  return (out)
}
k10 = function(del.l, del.u){
  out  = ((4*((del.l -1)**3))/((2*del.l -3)*(2*del.l -3)*(2*del.l -1)))
  return (out)
}
k11 = function(del.l, del.u){
  out  =  ((2)/((1 + 2*del.l)*(2*del.l -3)))
  return (out)
}
k12 = function(del.l, del.u){
  out  = (((-12*del.l + 12*(del.l**2) -5))/((1+2*del.l)*(1+2*del.l)*(2*del.l -3)*(2*del.l -3)))
  return (out)
}
##########marginal used to find inv#######################
pnormasympl <- function(x, del.l, del.u){
  if((del.l == 0.5) & (del.u == 0.5)){
    val <- 0.5 * (1 - x) * exp(2 * x) * as.numeric(x <= 0) + 
      (1 - 0.5 * (1 + x) * exp(-2 * x)) * as.numeric(x > 0)
  }else if((del.l == 0.5) & (del.u != 0.5)){
    val <- (k5(del.l, del.u) * x * exp(2 * x) - k6(del.l, del.u) * exp(2 * x)) * as.numeric(x <= 0) + 
      (1 - k7(del.l,del.u) * exp(-x / del.u) - k8(del.l, del.u) * exp(x / (del.u - 1))) * as.numeric(x > 0)
  }else if((del.l != 0.5) & (del.u == 0.5)){
    val <- (k9(del.l,del.u) * exp(x / del.l) + k10(del.l,del.u) * exp(x / (1 - del.l))) * as.numeric(x <= 0) + 
      (1 + k11(del.l,del.u) * x * exp(-2 * x) + k12(del.l, del.u) * exp(-2 * x)) * as.numeric(x > 0)
  }else{
    val <- (k1(del.l, del.u) * exp(x / del.l) - k2(del.l, del.u) * exp(x / (1 - del.l))) * as.numeric(x <= 0) + 
      (1 + k3(del.l, del.u) * exp(-x / del.u) - k4(del.l, del.u) * exp(x / (del.u-1))) * as.numeric(x > 0)
  }
  return(val)}

pnormasympl1 = function(x, del.l , del.u,t){
 solve = function(x, del.l , del.u,t){
   if(del.l == 0.5  && del.u == 0.5){
     if(x<=0){
       val =  abs(log(t) - log(0.5) - log((1-x)) - (2*x))
     }
     else {
       val = abs(log(1-t)  - log(0.5)  - log((1+x))  + (2*x))
    }
     return (val)
   }
#exp(1)
#else
   if(del.l == 0.5){
     if(x<=0){
       val =  abs(log(t)-   (2*x)  - log((k5(del.l,del.u)*x)  -  (k6(del.l,del.u))))
     }
    else {
       if(del.u > 0.5){
         val = abs(log(1-t)  + (x/del.u) - log(k7(del.l,del.u) + k8(del.l,del.u)*exp(x*((1-2*del.u)/(del.u*(1-del.u)))) ) )
       }
       else {
         val =  abs(log(1-t) + (x/(1-del.u)) - log(k7(del.l,del.u)*exp(x*((2*del.u -1)/(del.u*(1-del.u))))  + k8(del.l,del.u)) )
       }

     }
     return (val)
   }
   if(del.u == 0.5){
     if(x<=0){
       if(del.l >0.5){
         val = abs(log(t) - (x/del.l)  - log(k9(del.l,del.u) + k10(del.l,del.u)*exp(x*((2*del.l -1)/(del.l*(1-del.l))))) )
       }
       else {
         val = abs(log(t) - (x/(1-del.l))  - log(k9(del.l,del.u)*exp(x*((1-2*del.l)/(del.l*(1-del.l)))) + k10(del.l,del.u)) )
       }
     }
     else {
       val  = abs(log(1-t) + (2*x) - log(-1*(k11(del.l,del.u)*x + k12(del.l,del.u))) )
     }
     return (val)
   }
#else
   if(x<=0){
     if(del.l >0.5){
       val = abs(log(t) - (x/del.l) - log(k1(del.l,del.u) - k2(del.l,del.u)*exp(x*((2*del.l -1)/(del.l*(1-del.l)))) )  )
     }
     else {
       val = abs(log(t) - (x/(1-del.l)) - log(k1(del.l,del.u)*exp(x*((1-2*del.l)/(del.l*(1-del.l)))) - k2(del.l,del.u))  )
     }
   }
   else {
     if(del.u >0.5){
       val= abs(log(1-t) + (x/del.u) - log(-1*k3(del.l,del.u) + k4(del.l,del.u)*exp(x*((1-2*del.u)/(del.u*(1-del.u))))) )
     }
     else {
       val = abs(log(1-t) + (x/(1-del.u)) - log(-1*k3(del.l,del.u)*exp(x*((2*del.u-1)/(del.u*(1-del.u)))) + k4(del.l,del.u))  )
     }
 }
   return (val)
 }
 return (sapply(x,solve , del.l = del.l, del.u = del.u,t =t))
}



#x =seq(-4 , 100  ,by = 0.1)
#plot(x, pnormasympl(x))

#it is derivative of pnormasympl
#  f.dash  = function(x){


#      if(del.l ==0.5 && del.u ==0.5){
#      if(x<=0){
#        val =   -1* 0.5*(1-2*x)*exp(2*x)
#      }
#      else {
#       val =  -1*0.5*(1+2*x)*exp(-2*x)
#      }
#      return (val)
#      }

#       # else 
#    if(del.l == 0.5){
#        if(x<=0){ 
#          val =  t-   (exp(2*x)*2*x)/((1+2*del.u)*(2*del.u -3))  +  (exp(2*x)*(-12*del.u + 12*del.u*del.u  -5))/((1+2*del.u)*(1+2*del.u)*(2*del.u - 3)*(2*del.u - 3))
#        }
#        else {
#          val =  t -  1 +  (exp(-x/del.u)*4*(del.u**3))/((1+2*del.u)*(1+2*del.u)*(2*del.u -1)) + (exp(-x/(1-del.u))*4*((del.u -1)**3))/((del.u-3)*(del.u-3)*(2*del.u -1))
#        }
#        return (val)
#    }
#    #if(del.u ==.5){ }
#  # else 
#     if(x<=0){
#        val =  -1*((k1(del.l,del.u)*exp(x/del.l))/del.l) + ((k2(del.l,del.u)*exp(x/(1-del.l)))/(1-del.l))
#     }
#     else {
#        val =  ((k3(del.l,del.u)*exp(-x/del.u))/del.u) -  ((k4(del.l,del.u)*exp(x/(del.u-1)))/(1-del.u))
#     }
#     return (val)
#  }
#logit
logit= function(x, del.l,del.u,t){
  val =  abs(t- pnormasympl(log(x) - log(1-x),del.l,del.u))
  return (val)
}

#use this or newton raphson algorithm
# was giving error
#qasympl.sum = function(t,del.l ,del.u ){
# solve = function(t,del.l,del.u){
#    x = optimize(function(x, del.l , del.u,t){exp(pnormasympl(x, del.l , del.u,t))/(1+exp(pnormasympl(x, del.l , del.u,t)))} , c(0,1),tol = 1e-3, del.l = del.l, del.u = del.u, t = t)$minimum
#    return (x)
#  }
#  return (sapply(t, solve, del.l = del.l ,del.u = del.u))
#}
############
#qasympl.sum  with logit transformation applied
qasympl.sum = function(t,del.l ,del.u){
  solve = function(t,del.l,del.u){
    x = optimize(logit , c(0,1),tol = 1e-6, del.l = del.l, del.u = del.u, t = t)$minimum
    return (log(x) - log(1-x))
  }
  return (sapply(t, solve, del.l = del.l ,del.u = del.u))
}
#qasympl.sum(0.5, 0.2,0.3)
#but it is giving wrong ans for interval (1e-6, 1e6) what to do????????????????
#newton raphson
# tol = 1000
# epsilon = 1e-4
# iter =0
# root = 13
# while(tol > epsilon && iter <1000){
#     iter = iter +1
#     curr = root - pnormasympl(root , 0.2, 0.3)/f.dash(root)
#     tol = abs(root - curr)
#     root = curr
# }
# root
#integrand 
inte = function(x,x1 ,x2,rho,del.l,del.u){ 
  # passing log(x) -log(1-x) for logit transformation in integration to change limit of integration to (0,1 )
  solve = function(x ,x1 ,x2,rho,del.l,del.u){
    r <- log(x)-log(1-x)
    d1 <- dasymlp(x1 - r , 1 - del.l, 1 - del.u)
    q1 <- ifelse(d1 < 0.5, qnorm(d1), -qnorm(1 - d1))
    
    d2 <- dasymlp(x2 - r , 1 - del.l, 1 - del.u)
    q2 <- ifelse(d2 < 0.5, qnorm(d2), -qnorm(1 - d2))
    
    #q = c(qnorm(dasymlp(x1 - r , 1- del.l, 1-del.u)) , qnorm(dasymlp(x2-r , 1- del.l, 1-del.u)))
    q <- c(q1, q2)
    gcop =  pmvnorm(mean = c(0,0),corr  = matrix(c(1,rho,rho,1), nrow =2),lower = -Inf,upper = q)
    igrand = gcop[1]*denasymlp(r,del.l ,del.u)
    igrand = igrand/(x*(1-x)) # this comes due to substitution of r using logit transformation
    return (igrand)
  }
  return(sapply(x,solve,x1 =x1 ,x2 =x2, rho = rho, del.l =del.l ,del.u =del.u))
}

#pmvnorm(mean = c(0,0),corr  = matrix(c(1,0.3,0.3,1), nrow =2),lower = -Inf,upper = c(1,2))

# integration
C.t.t  = function(t1,t2,rho,del.l,del.u) {
  x1 = qasympl.sum(t1,del.l, del.u)
  x2 = qasympl.sum(t2,del.l , del.u)
  #Integration gives  Joint Distribution F(x1,x2)
  solve = function(x1,x2,rho,del.l,del.u){
    return (integrate(inte, lower = 0, upper = 1, x1 =x1 ,x2 =x2, rho = rho,del.l = del.l ,del.u =del.u)$value)
  }
  return (mapply(solve, x1,x2, rho =rho , del.l= del.l ,del.u =del.u))
}
#C.t.t(0.1,0.1,0.5,0.3,0.3)


#likelihood to estimate parameters
#1. empirical distribution to calc (u1,u2)
#x1 = c(3,4,3,5,3,4,2)
#x2 = c(3,5,5,6,7,5,4)
#pdat1 <- function(x,x1){
#  return(sum(x1 <= x)/(length(x1) +1))
#}
#pdat2 <- function(x,x2){
#  return(sum(x2 <= x)/(length(x2) +1))
#}
#2. c(u1,u2) 

#-----------------------

del.l = 0.3
del.u = 0.3
rho = 0.5

qasymlp <- function(p, del.l, del.u){
  cutoff <- del.l / (del.l + del.u)
  part1 <- (del.l + del.u) / (del.l * as.numeric(p <= cutoff) + del.u * as.numeric(p > cutoff))
  part2 <- ifelse(p <= cutoff, log(p), log(1 - p))
  part3 <- (del.l * as.numeric(p <= cutoff) - del.u * as.numeric(p > cutoff))
  out <- part3 * (log(part1) + part2)
  out}

#---------------

R <- qasymlp(runif(200), del.l, del.u)

Z <- rmvnorm(200, mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))

U <- pnorm(Z)

W <- qasymlp(U, 1 - del.l, 1 - del.u)

X <- W + R

 X <- apply(X, 2, rank) / (dim(X)[1]+1)
#X <- cbind(pnormasympl(X[,1],del.l,del.u),pnormasympl(X[,2],del.l,del.u)) 

#-----------------

# q <- cbind(qnorm(dasymlp(x1 - r , 1 - del.l, 1 - del.u)),
#            qnorm(dasymlp(x2 - r , 1 - del.l, 1 - del.u)))

den.integrand <- function(x, x1, x2, rho, del.l, del.u){ 
  r <-  log(x) - log(1 - x)  #logit transformation
  d1 <- dasymlp(x1 - r , 1 - del.l, 1 - del.u)
  q1 <- ifelse(d1 < 0.5, qnorm(d1), -qnorm(1 - d1))
  
  d2 <- dasymlp(x2 - r , 1 - del.l, 1 - del.u)
  q2 <- ifelse(d2 < 0.5, qnorm(d2), -qnorm(1 - d2))
  
  q <- cbind(q1, q2)
  
  out <- dmvnorm(q, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  out <- out * denasymlp(x1 - r, 1 - del.l, 1 - del.u) *
    denasymlp(x2 - r, 1 - del.l, 1 - del.u) * denasymlp(r, del.l, del.u)
  out <- out / (dnorm(dasymlp(x1 - r, 1 - del.l, 1 - del.u)) *
                  dnorm(dasymlp(x2 - r, 1- del.l, 1 - del.u))) 
  out <- out/(x*(1-x))
  out}

d.joint <- function(x1, x2, rho, del.l, del.u){
  integrate(den.integrand, lower = 0, upper = 1,
            x1 = x1, x2 = x2, rho = rho, del.l = del.l, del.u = del.u)$value
}

d.marginal = function(x,del.l,del.u){
  solve = function(x,del.l,del.u){
    #more cases of del.l = 1/2 and del.u =1/2 and .... 
    if(del.l == 0.5 && del.u == 0.5){
      if(x<=0){
        out =   ((1-x)*exp(2*x)) - (0.5*exp(2*x))
        return (out)
      } 
      else {
        out = ((1+x)*exp(-2*x)) - (0.5*exp(-2*x))
        return (out)
      }
    }
    # complete it
    if(del.l ==0.5){
      if(x<=0){
        out = (k5(del.l ,del.u)*(exp(2*x)*(1+2*x))) - (k6(del.l,del.u)*2*exp(2*x))
        return (out)
      }
      else {
        out = (((k7(del.l ,del.u)/del.u)*exp(-x/del.u))) + (((k8(del.l,del.u)/(1-del.u))*exp(-x/(1-del.u))))
        return (out)
      }
    } 
    if(del.u ==0.5){
      if(x<=0){
        out  = (((k9(del.l,del.u)/del.l)*exp(x/del.l))) + (((k10(del.l ,del.u)/(1-del.l))*exp(x/(1-del.l))))
        return (out)
      }
      else {
        out = (k11(del.l, del.u)*exp(-2*x)*(1-2*x)) - ((2*k12(del.l, del.u)*exp(-2*x)) )
        return (out)
      }
    }
    ##########
    if(x<=0){ 
      out = ( (k1(del.l,del.u)/del.l)*exp(x/del.l)) - ((k2(del.l,del.u)/(1-del.l))*exp(x/(1-del.l)) )
      return (out)
    }
    else {
      out = ((-k3(del.l,del.u)/del.u)*exp(-x/del.u)) +   ((k4(del.l,del.u)/(1-del.u))*exp(x/(del.u -1)))
      return (out)
    }
  }
  return (sapply(x,solve,del.l =del.l, del.u =del.u))
}

c.u1.u2 <- function(u1, u2, rho, del.l, del.u){
  x.1 <- qasympl.sum(u1, del.l, del.u)
  x.2 <- qasympl.sum(u2, del.l, del.u)
  out <- d.joint(x.1, x.2,rho,del.l , del.u) / 
    (d.marginal(x.1, del.l, del.u) * d.marginal(x.2, del.l, del.u))
  out}

#3. likelihood
neg.log <- function(params){
  sum(sapply(1:nrow(X), function(i){
    - log(c.u1.u2(X[i, 1], X[i, 2], #pdat1(x1[i],x1),pdat2(x2[i],x2),
                  (exp(params[1])-1)/(exp(params[1]) +1),
                  1/(1+exp(-params[2])),
                  1/(1+exp(-params[3])) ))}))}

 out <- optim(c(1.098612,-0.8472979,-0.8472979), neg.log,  hessian = TRUE, control = list(maxit = 1000))
 
 z <- NULL 
 
 z$nllh <- out$value
 z$conv <- out$convergence
 z$mle <- out$par
 z$cov <- solve(out$hessian)
 z$se <- sqrt(diag(z$cov))

# optim(c(0,0,0), neg.log)



#estimating lower tail cofficient
t = seq(0.01,0.999, by = 0.01)
x.l = numeric(length(t))
for(i in 1:length(t)){
  x.l[i] = C.t.t(t[i],t[i],rho,del.l ,del.u)/t[i]
} 
plot(t,x.l, type ='l')

#estimating upper tail coefficent
x.u = numeric(length(t))
for(i in 1:length(t)){
  x.u[i] = (1 - (2*t[i]) + C.t.t(t[i],t[i],rho,del.l,del.u))/(1 - t[i]);
} 
plot(t, x.u, type = 'l')

#plotting C.t.t 
install.packages("ggplot2")
library(ggplot2)
library(viridis)
loc = expand.grid(x=t, y=t)

fill = mapply( C.t.t,loc[1], loc[2], rho = 0.5, del.l = 0.1 , del.u = 0.9)
p1 = ggplot() + geom_tile(aes(x = loc[ , 1], y = loc[ , 2], fill = fill),
                          width = 0.06, height = 0.06) +
  coord_fixed(ratio = 1) +  scale_fill_viridis()


fill = mapply( C.t.t,loc[1], loc[2], rho = 0.5, del.l = 0.1 , del.u = 0.1)
p2 = ggplot() + geom_tile(aes(x = loc[ , 1], y = loc[ , 2], fill = fill),
                          width = 0.06, height = 0.06) +
  coord_fixed(ratio = 1) +  scale_fill_viridis()

fill = mapply( C.t.t,loc[1], loc[2], rho = 0.5, del.l = 0.5 , del.u = 0.5)
p3 = ggplot() + geom_tile(aes(x = loc[ , 1], y = loc[ , 2], fill = fill),
                          width = 0.06, height = 0.06) +
  coord_fixed(ratio = 1) +  scale_fill_viridis()

fill = mapply( C.t.t,loc[1], loc[2], rho = 0.5, del.l = 0.9 , del.u = 0.9)
p4 = ggplot() + geom_tile(aes(x = loc[ , 1], y = loc[ , 2], fill = fill),
                          width = 0.06, height = 0.06) +
  coord_fixed(ratio = 1)  +  scale_fill_viridis()


