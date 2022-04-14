
Indigo <- read_csv("INGL Historical Data.csv")
SpiceJet <- read_csv("SPJT Historical Data.csv")
Indigo <- Indigo[-c(497:507), ]
SpiceJet <- SpiceJet[-c(498:508), ]

years.Indigo <- as.vector(sapply(Indigo$Date, function(x){
  as.numeric(strsplit(x, split = " ")[[1]][3])}))

years.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){
  as.numeric(strsplit(x, split = " ")[[1]][3])}))

days.Indigo <- as.vector(sapply(Indigo$Date, function(x){
  as.numeric(strsplit(strsplit(x, split = ", ")[[1]][1], split = " ")[[1]][2])}))

days.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){
  as.numeric(strsplit(strsplit(x, split = ", ")[[1]][1], split = " ")[[1]][2])}))

months.Indigo <- as.vector(sapply(Indigo$Date, function(x){strsplit(x, split = " ")[[1]][1]}))

months.SpiceJet <- as.vector(sapply(SpiceJet$Date, function(x){strsplit(x, split = " ")[[1]][1]}))

months.Indigo <- as.vector(sapply(months.Indigo, function(x){which(month.abb == x)}))

months.SpiceJet <- as.vector(sapply(months.SpiceJet, function(x){which(month.abb == x)}))

identify.issue <- which(sapply(1:507, function(i){
  sum((months.Indigo - months.SpiceJet[-i])^2) +
    sum((years.Indigo - years.SpiceJet[-i])^2) +
    sum((days.Indigo - days.SpiceJet[-i])^2)}) == 0)

SpiceJet <- SpiceJet[-identify.issue, ]
years <- years.Indigo

months <- months.Indigo

days <- days.Indigo

nt <- length(Indigo$Date)

logreturn.Indigo <- log(Indigo$Price[-1]) - log(Indigo$Price[-nt])
logreturn.SpiceJet <- log(SpiceJet$Price[-1]) - log(SpiceJet$Price[-nt])

par(mfrow = c(1, 2))
acf(logreturn.Indigo)
acf(logreturn.SpiceJet)

######################################
dates <- mdy(Indigo$Date)

df <- data.frame(date = dates[-1], 
                 Indigo =logreturn.Indigo,
                 SpiceJet =logreturn.SpiceJet)

##############################################
p <- ggplot(df) + geom_line(aes(x = date, y = SpiceJet), size = 1) + 
  xlab(NULL) + ylab("Average-Price") + 
  ggtitle("Daily average price of SpiceJet airlines") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))

print(p)
###########################################
p <- ggplot(df) + geom_point(aes(x = Indigo, y = SpiceJet)) + ggtitle("Bivariate Plot of SpiceJet and Indigo logreturn data") +
  theme(plot.title = element_text(hjust = 1))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))
print(p)

################################################
Indigo.rank <- (rank(df$Indigo) + 0.5) / (length(df$Indigo) + 1)
SpiceJet.rank <- (rank(df$SpiceJet) + 0.5) / (length(df$SpiceJet) + 1)

df2 <- data.frame(Indigo.rank = Indigo.rank,
                  SpiceJet.rank = SpiceJet.rank)

p <- ggplot(df2) + geom_point(aes(x = Indigo.rank, y = SpiceJet.rank)) + ggtitle("INGL and SPJT scatter plot after rank transformation") +
  theme(plot.title = element_text(hjust = 1))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=20))
print(p)

##################################################

C.t.t = function(t,del.l,del.u,rho){
  solver <- function(t,del.l,del.u,rho){
    x<- qasympl.sum(t,del.l,del.u)
    out <- pasympl.gauss1(t,t, del.l,del.u,rho)
    out
  }
  sapply(t, solver , del.l = del.l ,del.u = del.u,rho=rho)
}
chi.l <- function(t, del.l,del.u,rho){
  out = C.t.t(t,del.l,del.u,rho)/t
  out
} 

chi.u <- function(t,del.l ,del.u,rho){
  out <- (1 - (2*t) + C.t.t(t,del.l,del.u,rho))/(1-t)
  out
}
####################################### 
rho = mles[3]
del.l = mles[1] 
del.u = mles[2]

base <- ggplot() + xlim(0.1, 0.9)

base + geom_function(fun = chi.u, args = list(del.l  = del.l, del.u = del.u,rho = rho))+ 
               # coord_cartesian(ylim = c(0, 1)) +
  xlab("x") + ylab("Lower tail Coefficient")
