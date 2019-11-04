# ============================= #
#        News Factor
# Author : Kyounghoon Kim
# Date : 2019.09.19
# ============================= #

# Change

# For graphics

install.packages("graphics")
require(graphics)

# Test Cross_correlation
install.packages("corrplot")
library(corrplot)

# Stationary test
install.packages("xts")
install.packages("tseries")
install.packages("Quandl")
install.packages("urca")
library('xts')
library('tseries')
library('Quandl')
library('urca')

# For Time Series
install.packages("TTR")
install.packages("forecast")
install.packages("stringe")
library(TTR)
library(forecast)
library(stringr)

# For Distance
install.packages("dtw")
library(dtw)


# Functions
# Stationary test Looks have Seasonality
#

Season.input.func <- function(data.set,time.model) # auto.arima로 나온 결과를 활용해 seasonal 잡아주는 함수. seasonal order 1 or 2
{
  a = c(1,2)
  b <- matrix(c(1,1,combn(a,2),2,1,2,2),nrow=2)
  
  best.model <- time.model
  ar.order <- time.model$arma[1]
  ma.order <- time.model$arma[2]
  arima.diff <- time.model$arma[5]
  Period_list <- c(5,7,12)
  
  for (i in (1:4))
  {
    seasonal.ar <- b[,i][1]
    seasonal.ma <- b[,i][2]
    temp.model <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,0,seasonal.ma),period = 5),method="CSS")
    temp.model2 <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,0,seasonal.ma),period = 7),method="CSS")
    temp.model3 <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,0,seasonal.ma),period = 12),method="CSS")
    
    temp.model4 <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,1,seasonal.ma),period = 5),method="CSS")
    temp.model5 <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,1,seasonal.ma),period = 7),method="CSS")
    temp.model6 <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,1,seasonal.ma),period = 12),method="CSS")
    
    min.index <- which.min(c(temp.model$loglik*-2,temp.model2$loglik*-2,temp.model3$loglik*-2,temp.model4$loglik*-2,temp.model5$loglik*-2,temp.model6$loglik*-2))
    if( min.index <= 3 )
    {
      Season.diff <- 0 
      Can.Period <- Period_list[min.index]
      cand.model <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,Season.diff,seasonal.ma),period = Can.Period),method="CSS")
      if( best.model$loglik*-2 >= cand.model$loglik*-2)
      {
        best.model <- cand.model
      }
    }
    else
    {
      Season.diff <- 1
      Can.Period <- Period_list[min.index-3]
      cand.model <- arima(data.set,order=c(ar.order,arima.diff,ma.order),seasonal=list(order = c(seasonal.ar,Season.diff,seasonal.ma),period = Can.Period),method="CSS")
      if( best.model$loglik*-2 >= cand.model$loglik*-2)
      {
        best.model <- cand.model
      } 
    }
    print(i)
  }
  return(c(ar.order,arima.diff,ma.order,seasonal.ar,Season.diff,seasonal.ma,Can.Period))
}

Stationary.func <- function(data_set,col.name,print.plot = "no")
{
  set <- data_set
  name <- as.character(col.name)
  # ================================ #
  # 안정화 과정 #
  Raw.model <- auto.arima(set[,name])
  Serial <- Season.input.func(set[,name],Raw.model)
  result <- arima(set[,name],order=c(Serial[1],Serial[2],Serial[3]),seasonal = list(order=c(Serial[4],Serial[5],Serial[6]),period=Serial[7]),method="CSS")
  new1 <- result$residuals
  
  if ( print.plot == "yes")
  {
    # 안정화 전 시계열 acf, pacf 저장
    Before.s.path <- "C:/Users/code5/OneDrive/바탕 화면/KSIF(크시프)/뉴스경기판단/시계열안정화플롯과정/안정화전/"
    
    # acf 저장
    acf.path <- paste(Before.s.path,"/",name,".acf",".jpg",sep="")
    jpeg(acf.path)
    plot(acf(set[,name]))
    dev.off()
    
    # pacf 저장
    pacf.path <- paste(Before.s.path,"/",name,".pacf",".jpg",sep="")
    jpeg(pacf.path)
    plot(pacf(set[,name]))
    dev.off()
    
    # 잔차 플롯
    resid.path <- paste(Before.s.path,"/",name,".resid",".jpg",sep="") 
    jpeg(resid.path)
    plot(set[,name],type="l")
    dev.off()
    
    # 안정화 후 시계열 acf, pacf 저장
    after.s.path <- "C:/Users/code5/OneDrive/바탕 화면/KSIF(크시프)/뉴스경기판단/시계열안정화플롯과정/안정화후/"
    
    # acf 저장
    acf.path2 <- paste(after.s.path,"/",name,".acf",".jpg",sep="")
    jpeg(acf.path2)
    plot(acf(new1))
    dev.off()
    
    # pacf 저장
    pacf.path2 <- paste(after.s.path,"/",name,".pacf",".jpg",sep="")
    jpeg(pacf.path2)
    plot(pacf(new1))
    dev.off()
    
    # 잔차 플롯
    resid.path2 <- paste(after.s.path,"/",name,".resid",".jpg",sep="") 
    jpeg(resid.path2)
    plot(new1,type="l")
    dev.off()
    
    print("Export Plot")
    
  }
  
  return(as.numeric(new1))
}



# Sector Index
SECTOR_INDEX <- read.csv(paste(New.path,"SECTOR_INDEX2.csv",sep=""),header=T)

Index.name_list <- names(SECTOR_INDEX)[2:length(SECTOR_INDEX)]
SECTOR_residual <- data.frame(nrow=nrow(SECTOR_INDEX))

for ( col.name in Index.name_list)
{
  SECTOR_Resid <- Stationary.func(SECTOR_INDEX,col.name,"no") # yes: Export plot, no : return residual  
  SECTOR_residual <- cbind(SECTOR_residual,SECTOR_Resid)
  print(paste(col.name," Process Done",sep=""))
}
SECTOR_residual <- SECTOR_residual[,-1]
names(SECTOR_residual) <- Index.name_list   # Make Stationary process of data
SECTOR_residual


# ========================================================================================== #
#      News Data , Daily
#   Make Stationary For CCF
# ========================================================================================== #


path.news.data <- ("C:/Users/code5/OneDrive/바탕 화면/KSIF(크시프)/뉴스경기판단/빅카인즈 데이터/")
news.file.list <- list.files(path.news.data)
news.file.list[1] # 일별데이터
Daily.news.data.path <- paste(path.news.data,news.file.list[1],sep="")
Daily.news.file.list <- list.files(Daily.news.data.path)
Daily.news.file.list
News.data.set <- read.csv(paste(Daily.news.data.path,"/",Daily.news.file.list[1],sep=""),header=F)
News.data.set<-News.data.set[-1,1:2]
names(News.data.set) <- c("Date",str_sub(Daily.news.file.list[1],1,-5))


for( i in (2:length(Daily.news.file.list)))
{
  temp.set <- read.csv(paste(Daily.news.data.path,"/",Daily.news.file.list[i],sep=""),stringsAsFactors = F,encoding = "UTF-8")
  names(temp.set) <- c("Date",str_sub(Daily.news.file.list[i],1,-5))
  News.data.set <- merge(News.data.set,temp.set,by="Date")
}
News.data.set[,2] <- as.numeric(as.character(News.data.set[,2]))
News.data.set <- News.data.set[,-ncol(News.data.set)]
str(News.data.set)



News.name_list <- names(News.data.set)[2:length(News.data.set)]
News.name_list
Sta.news.data <- data.frame(nrow=nrow(News.data.set))


for ( col.name in News.name_list)
{
  Resid <- Stationary.func(News.data.set,col.name,"no") # yes: Export plot, no : return residual  
  Sta.news.data <- cbind(Sta.news.data,Resid)
  print(paste(col.name," Process Done",sep=""))
}
ncol(Sta.news.data)
Sta.news.data <- Sta.news.data[,-1]
names(Sta.news.data) <- News.name_list   # Make Stationary process of data


# ===================================================== 
#                  Find CCF
# =====================================================
# Compare with SECTOR_residual

install.packages("gtools")
library(gtools)

sector_col_num <- ncol(SECTOR_residual)
news_col_num <- ncol(Sta.news.data)

# Cross Correlation Interval
n <- nrow(Sta.news.data)
upper.bound <- -(1/n) + (2/sqrt(n))
lower.bound <- -(1/n) - (2/sqrt(n))


Relation_matrix <- matrix(0,nrow=sector_col_num,ncol = news_col_num)
Relation_matrix
rownames(Relation_matrix) <- Index.name_list
colnames(Relation_matrix) <- News.name_list

for ( i in (1:sector_col_num))
{
  for ( j in (1:news_col_num))
  {
    temp <- ccf(SECTOR_residual[,i],Sta.news.data[,j])
    temp.vector <- as.vector(temp$lag)
    
    # point out positive or negative
    if( (abs(max(temp$acf)) > abs(min(temp$acf))) && (max(temp$acf) >= upper.bound) ) # positive
    {
      Relation_matrix[i,j] <- temp.vector[temp$acf == max(temp$acf)]
    }
    else if ((abs(max(temp$acf)) < abs(min(temp$acf))) && (min(temp$acf) <= lower.bound) )
    {
      Relation_matrix[i,j] <- temp.vector[temp$acf == min(temp$acf)]
    }
    
  }
}

Relation_matrix

# Find 7 or -7 value of ccf

CCF_index_1 <- which(Relation_matrix == 7,arr.ind = T)
name.array_1 <- list()

for ( i in (1:nrow(CCF_index_1)))
{
  v <- CCF_index_1[i,]
  row_name <- rownames(Relation_matrix)[v[1]]
  col_name <- colnames(Relation_matrix)[v[2]]
  name.array_1[[i]] <- c(row_name,col_name)  
}
name.array_1

CCF_index_2 <- which(Relation_matrix == -7,arr.ind = T)
name.array_2 <- list()

for ( i in (1:nrow(CCF_index_2)))
{
  v <- CCF_index_2[i,]
  row_name <- rownames(Relation_matrix)[v[1]]
  col_name <- colnames(Relation_matrix)[v[2]]
  name.array_2[[i]] <- c(row_name,col_name)  
}
name.array_2

# =============================== #
# DTW clustering !!@!@!@!@!@!@!@
# 2019.9월 자로 새로 신선하게 나온 Package!
# =============================== #
install.packages("dtwclust")
library(dtwclust)
install.packages("ggdendro")
library(ggdendro)
library(dplyr)
library(ggplot2)
library(gtable)
library(gridExtra)



# vignette("dtwclust") dtw paper pdf
ndtw <- function(x, y, ...) {
  dtw(x, y, ...,
      step.pattern = asymmetric,
      distance.only = TRUE)$normalizedDistance  # Extract NormalizedDistance
}

# Register the distance with proxy
proxy::pr_DB$set_entry(FUN = ndtw, names = c("nDTW"),
                       loop = TRUE, type = "metric", distance = TRUE,
                       description = "Normalized, asymmetric DTW")

Dtw_clust <- dtwclust::tsclust(t(News.data.set[,2:ncol(News.data.set)]), 
                               type = "h", # hierarchical Clustering
                               k = 3,  # Number of Groups
                               distance = "nDTW", # distance method
                               control = hierarchical_control(method = "complete"),
                               preproc = NULL, 
                               args = tsclust_args(dist = list(window.size = 5L)))

hclus <- stats::cutree(Dtw_clust, k = 3) %>% # hclus <- cluster::pam(dist_ts, k = 2)$clustering has a similar result
  as.data.frame(.) %>%
  dplyr::rename(.,cluster_group = .) %>%
  tibble::rownames_to_column("type_col")

hcdata <- ggdendro::dendro_data(Dtw_clust)
names_order <- hcdata$labels$label
# Use the folloing to remove labels from dendogram so not doubling up - but good for checking hcdata$labels$label <- ""

p1 <- hcdata %>%
  ggdendro::ggdendrogram(., rotate=TRUE, leaf_labels=FALSE)

p2 <- Sta.news.data %>%
  dplyr::mutate(index = 1:987) %>%
  tidyr::gather(key = type_col,value = value, -index) %>%
  dplyr::full_join(., hclus, by = "type_col") %>% 
  mutate(type_col = factor(type_col, levels = rev(as.character(names_order)))) %>% 
  ggplot(aes(x = index, y = value, colour = cluster_group)) +
  geom_line() +
  facet_wrap(~type_col, ncol = 1, strip.position="left") + 
  guides(color=FALSE) +
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_blank())

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2) 

grid.arrange(gp2, gp1, ncol=2, widths=c(4,2))


















