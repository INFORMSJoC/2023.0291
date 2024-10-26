library(Hmisc)
data <- read.csv("data.csv", header=TRUE,
                             fileEncoding = "GBK",encoding = "GBK" )
names(data)
data<-data[,-c(1,2)]
data<-na.omit(data)
data<-data[-which(data$揽收量 %in% 0),]
y<-data[,1]
y<-log(y)
X<-data[,-1]
X$派件量<-log(X$派件量)

################# sort by correlation #################
cor<-cor(X,y)
X.ordered<-X[,order(abs(cor), decreasing = TRUE)]
data_sampled<-data.frame(y,X.ordered)
write.csv(data_sampled,'data_sampled.csv',fileEncoding = "GBK")

cor.ordered<-data.frame(Correlation=cor[order(abs(cor), decreasing = TRUE),])
cor.ordered$Names<-row.names(cor.ordered)
cor.ordered$Names<-factor(cor.ordered$Names,levels = cor.ordered$Names)
row.names(cor.ordered)<-1:dim(cor.ordered)[1]
write.csv(cor.ordered,'Correlation.csv',fileEncoding = "GBK")

################# sort by p-value #################
model<-lm(y~., data = X)
p_value<-summary(model)$coefficients[-1,4]
covariate<-X[,order(p_value, decreasing = FALSE)]
result<-data.frame(covariate=colnames(covariate), 
                   p_value=p_value[order(p_value, decreasing = FALSE)])
rownames(result)<-1:dim(result)[1]
write.csv(result, 'p_value.csv',fileEncoding = "GBK")
write.csv(data.frame(y,covariate), 'data_sampled2.csv',fileEncoding = "GBK")



