library(Hmisc)
library(readxl)
library(caTools)
library(dplyr)

PredictorData <- read_excel("PredictorData2022.xlsx")
PredictorData <- as.data.frame(lapply(PredictorData, as.numeric))

data<-PredictorData |>mutate(
d.p = log(D12) - log(Index), # Dividend Price ratio
d.y = log(D12) - log(lag(Index)), # Dividend yield
e.p = log(E12) - log(Index), # Earnings price ratio
d.e = log(D12) - log(E12), # Dividend payout ratio
tms = lty - tbl, # Term spread
dfy = BAA - AAA, # Default yield spread
dfr = corpr - ltr, # Default Return Spread
y = lead(CRSP_SPvw-Rfree), # excess stock returns
se.p = runmean(E12,k=10)/Index, # Smoothed Earnings Price Ratio (se/p)
lagy = lag(y,1) # Lagged dependent variable
) |> select(yyyymm, y, d.y, d.p, e.p, b.m, ntis, tbl, lty, ltr, dfy, dfr, infl, # 2008
              se.p, svar, lagy) |> filter(yyyymm >= 195001 & yyyymm <= 201012) 
data <- na.omit(data)

y<-data$y
X<-data[,-which(names(data) %in% c('yyyymm','y') )]

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



