library (tidyverse)
library (caret)

#Download zipped csv to data
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/00372/HTRU2.zip'
file = 'HTRU_2.csv'
dl = tempfile()
download.file(url, dl)
con = unz(dl, file)
con = unz('HTRU2.zip', file)
data = read.csv(con, header = FALSE)
colnames(data) = c('pmean', 'psd', 'pkurt', 'pskew', 'cmean', 'csd', 'ckurt', 'cskew', 'pulsar')

#Change pulsar to factor
data$pulsar=as.factor(data$pulsar)

#Keep .1 of data for validation
ind = createDataPartition(data$pulsar, times = 1, p = .7, list = FALSE)
train = data[ind,]
test = data[-ind,]

#make bootstrap indicies
strapcount = 100
bootind = createResample(train$pulsar, times = strapcount, list = FALSE)

#Returns hit ratio and f score for given model
f_acc = function(model){
  f = 0
  hr = 0
  for (i in c(1:strapcount)){
    bootstrap = train[bootind[,i],]
    temp = bootstrap %>% mutate(pulsarhat = predict(model)) %>% 
      select(pulsar, pulsarhat) %>% mutate(hit = pulsar == pulsarhat)
    f = f + F_meas(temp[,'pulsarhat'], reference = temp[,'pulsar'], beta = 2)
    hr = hr + mean(temp$hit)
  }
  f = f/strapcount
  hr = hr/strapcount
  data.frame(f, hr)
}

#exploratory boxplots
explore = train %>% gather('key', 'value', 1:8)
explore %>%  ggplot(aes(x = pulsar, y = value, group = pulsar)) + 
  geom_boxplot()  + facet_wrap(~key, scales = 'free')

#Start with csd and weed out ones that have strong dependence
machinefood = train
explore = machinefood[,names(machinefood) != 'pulsar'] %>% gather('key', 'value', -one_of('csd'))
explore %>% 
  ggplot(aes(x = csd, y = value, group = key)) + geom_point() + facet_wrap(~key, scales = 'free')

machinefood = machinefood[, !names(train) %in% c('ckurt', 'cskew')]

#Add pkurt and weed out ones that have strong dependence
explore = machinefood[,names(machinefood) != 'pulsar'] %>% gather('key', 'value', -one_of('pkurt'))
explore %>% 
  ggplot(aes(x = pkurt, y = value, group = key)) + geom_point() + facet_wrap(~key, scales = 'free')

machinefood = machinefood[, !names(machinefood) %in% c('pmean', 'pskew')]

##########################Probably garbage, remanent from old strategy
# explore = train[,c(1,2,3,5,6,7,8)] %>% gather('key', 'value', c(1,2,4,5,6,7))
# explore %>%  ggplot(aes(x = pkurt, y = value, group = pkurt)) + 
#   geom_point()  + facet_wrap(~key, scales = 'free')
# 
# explore = train[,c(1,2,4,5,6,7,8)] %>% gather('key', 'value', c(1,2,4,5,6,7))
# explore %>%  ggplot(aes(x = pskew, y = value, group = pskew)) + 
#   geom_point()  + facet_wrap(~key, scales = 'free')
# 
# train %>% ggplot(aes(x = pkurt, y = pmean)) + geom_point()



#pkurt and pskew are good candidates
#because tight 50% on edge of range
#slope is almost zero 
#at low pkurt (which is where 1 & 0 overlap) and where most of the data is
#csd has good separation too and is not strongly correlated to kurt or skew
# train %>% ggplot(aes(x = pkurt, y = pskew)) +geom_point()
# train %>% ggplot(aes(x = csd, y = pskew)) +geom_point()
# train %>% ggplot(aes(x = pkurt, y = csd)) +geom_point()
# 
# train %>% select(pulsar, pkurt, pskew, csd) %>% 
#   gather('key', 'value', c(pkurt, pskew, csd)) %>%
#   ggplot(aes(x=pulsar, y = value, group = pulsar)) + 
#   geom_boxplot() + facet_wrap(~key, scales = 'free')
########################################################

#random forest model
rfmodel1 = train(pulsar ~ pkurt + pskew, train, method = 'Rborist')
rfmodel1stats = f_acc(rfmodel1)

rfmodel2 = train(pulsar ~ pkurt + pskew + csd, train, method = 'Rborist')
start = timestamp()
rfmodel2stats = f_acc(rfmodel2)
end = timestamp()

start = timestamp()
rfmodel3 = train(pulsar ~. , machinefood, method = 'Rborist')
rfmodel3stats = f_acc(rfmodel2)
end = timestamp()

saveRDS(rfmodel2, file = 'rfmodel2.Rda')
saveRDS(rfmodel1, file = 'rfmodel1.Rda')
saveRDS(rfmodel3, file = 'rfmodel3.Rda')
saveRDS(rfmodel1stats, file = 'rfmodel1stats.Rda')
saveRDS(rfmodel3stats, file = 'rfmodel3stats.Rda')
#


