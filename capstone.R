library (tidyverse)
library (caret)
set.seed(5)
#Download zipped csv to data
url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/00372/HTRU2.zip'
file = 'HTRU_2.csv'
dl = tempfile()
download.file(url, dl)
con = unz(dl, file)
con = unz('HTRU2.zip', file)
data = read.csv(con, header = FALSE)
colnames(data) = c('pmean', 'psd', 'pkurt', 'pskew', 'cmean', 'csd', 'ckurt', 'cskew', 'pulsar')

Beta = 2

#Change pulsar to factor
data$pulsar=as.factor(data$pulsar)

#Keep .3 of data for validation
ind = createDataPartition(data$pulsar, times = 1, p = .7, list = FALSE)
train = data[ind,]
test = data[-ind,]

#make bootstrap indicies
strapcount = 100
bootind = createResample(train$pulsar, times = strapcount, list = FALSE)

#Returns hit ratio and f score for given model
modelstats = function(model, b = Beta){
  f = 0
  sens = 0
  spec = 0
  for (i in c(1:strapcount)){
    bootstrap = train[bootind[,i],]
    temp = bootstrap %>% mutate(pulsarhat = predict(model, newdata = bootstrap)) %>%
      select(pulsar, pulsarhat)
    cm = confusionMatrix(temp$pulsarhat, reference = temp$pulsar)
    sens = sens + cm$byClass['Sensitivity']
    spec = spec + cm$byClass['Specificity']
    f = f + F_meas(temp[,'pulsarhat'], reference = temp[,'pulsar'], beta = b)
  }
  df = data.frame(c(sens, spec, f, b)) %>% t()
  colnames(df) = c('sensitivity', 'specificity','fScore', 'beta')
  df
}

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F_meas(data$pred, reference = data$obs, beta = Beta)
  c(F1 = f1_val)
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

#check if cmean and psd are not dependent
explore = machinefood[,c('cmean', 'psd')]
explore %>% ggplot(aes(x = cmean, y = psd)) + geom_point() 

rfmodel = train(pulsar~., 
                 machinefood, 
                 method = 'Rborist', 
                 trControl = trainControl(method = 'boot', summaryFunction = f1))
allstats = modelstats(rfmodel, Beta) %>% data.frame() %>% mutate(model = 'Rborist')

for (b in c(2:5)){
  Beta = b
  print(Beta)
  rfmodel = train(pulsar~., 
                  machinefood, 
                  method = 'Rborist', 
                  trControl = trainControl(method = 'boot', summaryFunction = f1))
  temp = modelstats(rfmodel, Beta) %>% data.frame() %>% mutate(model = 'Rborist')
  allstats = union(allstats,temp)
  
  knnmodel = train(pulsar~., 
                   machinefood, 
                   method = 'knn', 
                   trControl = trainControl(method = 'boot', summaryFunction = f1))
  temp = modelstats(knnmodel, Beta) %>% data.frame() %>% mutate(model = 'knn')
  allstats = union(allstats, temp)
  
  glmmodel = train(pulsar~., 
                   machinefood, 
                   method = 'glm', 
                   trControl = trainControl(method = 'boot', summaryFunction = f1))
  temp = modelstats(glmmodel, Beta) %>% data.frame() %>% mutate(model = 'glm')
  allstats = union(allstats, temp)
  
  
  nbmodel = train(pulsar~., 
                  machinefood, 
                  method = 'nb', 
                  trControl = trainControl(method = 'boot', summaryFunction = f1))
  temp = modelstats(nbmodel, Beta) %>% data.frame() %>% mutate(model = 'nb')
  allstats = union(allstats, temp)
  
  svmmodel = train(pulsar~., 
                   machinefood, 
                   method = 'svmLinearWeights2', 
                   trControl = trainControl(method = 'boot', summaryFunction = f1))
  temp = modelstats(svmmodel, Beta) %>% data.frame() %>% mutate(model = 'svm')
  allstats = union(allstats, temp)
}
rownames(stats) = c(1:100)


# saveRDS(machinefood, file = 'machinefood.Rda')
# saveRDS(rfmodel, file = 'rfmodel.Rda')
# saveRDS(rfmodel1stats, file = 'rfmodel1stats.Rda')
# saveRDS(rfmodel3stats, file = 'rfmodel3stats.Rda')
# saveRDS(knnmodel, file = 'knnmodel.Rda')
# saveRDS(knnmodelstats, file = 'knnmodelstats.Rda')
# saveRDS(glmmodel, file = 'glmmodel.Rda')
# saveRDS(glmmodelstats, file = 'glmmodelstats.Rda')
# saveRDS(nbmodel, file = 'nbmodel.Rda')
# saveRDS(nbmodelstats, file = 'nbmodelstats.Rda')
# saveRDS(svmmodel, file = 'svmmodel.Rda')
# rfmodel3 = read_rds('rfmodel3.Rda')
# knnmodel = read_rds('knnmodel.Rda')

