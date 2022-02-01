data1<- read.csv("ingl.csv")
data2<- read.csv("spicejet.csv")


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
data1$Date = as.Date(data1$Date, format="%b %d, %Y")
data1<-data1[order(as.Date(data1$Date, format="%b %d, %Y")),]
plot(data1$Date, data1$Price, type = 'l')

log1 = numeric(length = dim(data1)[1] -1)
for(i in 2:dim(data1)[1]){
  log1[i-1]  = log(data1$Price[i]) - log(data1$Price[i-1])
  
}
plot(log1, type = 'l')


##############data2 plots#################

data2$Date = as.Date(data2$Date, format="%b %d, %Y")
data2<-data2[order(as.Date(data2$Date, format="%b %d, %Y")),]
plot(data2$Date, data2$Price, type = 'l')
# log return
log2 = numeric(length = dim(data2)[2] -1)
for(i in 2:dim(data2)[2]){
  log2[i-1]  = log(data2$Price[i]) - log(data2$Price[i-1])
  
}
plot(log2, type = 'l')



