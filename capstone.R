library (tidyverse)
library (caret)

#Download zipped csv to data
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/00372/HTRU2.zip'
file = 'HTRU_2.csv'
dl = tempfile()
download.file(url, dl)
con = unz(dl, file)
data = read.csv(con, header = FALSE)
colnames(data) = c('pmean', 'psd', 'pkurt', 'pskew', 'cmean', 'csd', 'ckurt', 'cskew', 'pulsar')

#Change pulsar to factor
data$pulsar=as.factor(data$pulsar)

#Keep .1 of data for validation
ind = createDataPartition(data$pulsar, times = 1, p = .8, list = FALSE)
train = data[ind,]
test = data[-ind,]

#make bootstrap indicies
strapcount = 100
bootind = createResample(train$pulsar, times = strapcount, list = FALSE)

#exploratory boxplots
explore = train %>% gather('key', 'value', 1:8)
explore %>%  ggplot(aes(x = pulsar, y = value, group = pulsar)) + 
  geom_boxplot()  + facet_wrap(~key, scales = 'free')

#pkurt and pskew are good candidates
#because tight 50% on edge of range
#and pkurt and pskew are not strongly correlated
#at low pkurt (which is where 1 & 0 overlap)
train %>% ggplot(aes(x = pkurt, y = pskew)) +geom_point()
train %>% ggplot(aes(x=pulsar, y = pkurt, group = pulsar)) + geom_boxplot()
train %>% ggplot(aes(x=pulsar, y = pskew, group = pulsar)) + geom_boxplot()

#random forest model
rfmodel = train(pulsar ~ pkurt + pskew, train, method = 'Rborist')
saveRDS(rfmodel, file = 'rfmodel.Rda')

#Estimate of hit ratio
sum = 0
for (i in c(1:strapcount)){
  bootstrap = train[bootind[,i],]
  temp = bootstrap %>% mutate(pulsarhat = predict(rfmodel)) %>% 
    mutate(hit = pulsar == pulsarhat) %>%
    summarise(mean(hit))
  sum = sum + temp
}
mean = sum/strapcount



#


