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

#Change pulsar to factor
data$pulsar=as.factor(data$pulsar)

#Keep .3 of data for validation
ind = createDataPartition(data$pulsar, times = 1, p = .7, list = FALSE)
train = data[ind,]
test = data[-ind,]
saveRDS(train, 'train.rda')
saveRDS(test, 'test.rda')

#Make bootstrap indicies
strapcount = 100
bootind = createResample(train$pulsar, times = strapcount, list = FALSE)

#Returns fscore, sensitivity, and specificity for a given model by bootstrapping
modelstats = function(model, b = Beta){
  f = 0
  sens = 0
  spec = 0
  for (i in c(1:strapcount)){
    #Create bootstrap
    bootstrap = train[bootind[,i],]
    #Estimate outcome
    temp = bootstrap %>% mutate(pulsarhat = predict(model, newdata = bootstrap)) %>%
      select(pulsar, pulsarhat)
    cm = confusionMatrix(temp$pulsarhat, reference = temp$pulsar)
    #Add values to be averaged
    sens = sens + cm$byClass['Sensitivity']
    spec = spec + cm$byClass['Specificity']
    f = f + F_meas(temp[,'pulsarhat'], reference = temp[,'pulsar'], beta = b)
  }
  #Take averages
  sens = sens/strapcount
  spec = spec/strapcount
  f = f/strapcount
  #Put in data frame and rename
  df = data.frame(c(sens, spec, f, b)) %>% t()
  colnames(df) = c('sensitivity', 'specificity','fScore', 'beta')
  df
}
#Set Beta for F-score
Beta = 1

#Function to return F1 score, used for training models
f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F_meas(data$pred, reference = data$obs, beta = Beta)
  c(F1 = f1_val)
}

#Function to return specificity, used for training models
spec <- function(data, lev = NULL, model = NULL) {
  c(Specificity = specificity(data = data$pred, reference = data$obs))
}

#Function to return Sensitivity, used for training models
sens <- function(data, lev = NULL, model = NULL) {
  c(Sensitivity = sensitivity(data = data$pred, reference = data$obs))
}

#Exploratory boxplots
train %>% gather('key', 'value', 1:8) %>%  
  ggplot(aes(x = pulsar, y = value, group = pulsar)) + 
  geom_boxplot()  + facet_wrap(~key, scales = 'free')

#Start with csd and weed out ones that have strong dependence
train[,names(train) != 'pulsar'] %>% 
  gather('key', 'value', -one_of('csd')) %>% 
  ggplot(aes(x = csd, y = value, group = key)) + 
  geom_point() + facet_wrap(~key, scales = 'free')

machinefood = train
machinefood = machinefood[, !names(train) %in% c('ckurt', 'cskew')]

#Add pkurt and weed out ones that have strong dependence
train[,!names(train) %in% c('pulsar', 'ckurt', 'cskew')] %>% 
  gather('key', 'value', -one_of('pkurt'))%>% 
  ggplot(aes(x = pkurt, y = value, group = key)) + 
  geom_point() + facet_wrap(~key, scales = 'free')

machinefood = machinefood[, !names(machinefood) %in% c('pmean', 'pskew')]

#check if cmean and psd are dependent
train[,c('cmean', 'psd')]%>% ggplot(aes(x = cmean, y = psd)) + geom_point(size = .1) +
  ggtitle('Figure 4')
#Train models with default metric
rfmodel = train(pulsar~., 
                 machinefood, 
                 method = 'Rborist', 
                 trControl = trainControl(method = 'boot'))

knnmodel = train(pulsar~., 
                machinefood, 
                method = 'knn', 
                trControl = trainControl(method = 'boot'))

glmmodel = train(pulsar~., 
                 machinefood, 
                 method = 'glm', 
                 trControl = trainControl(method = 'boot'))

nbmodel = train(pulsar~., 
                 machinefood, 
                 method = 'nb', 
                 trControl = trainControl(method = 'boot'))

#Estimate performance of each model and add them to allstats data frame
allstats = modelstats(rfmodel, Beta) %>% data.frame() %>% mutate(model = 'Rborist no metric')

temp = modelstats(knnmodel, Beta) %>% data.frame() %>% mutate(model = 'knn no metric')
allstats = union(allstats,temp)

temp = modelstats(glmmodel, Beta) %>% data.frame() %>% mutate(model = 'glm no metric')
allstats = union(allstats,temp)

temp = modelstats(nbmodel, Beta) %>% data.frame() %>% mutate(model = 'nb no metric')
allstats = union(allstats,temp)

#Train models with varying F1 metric and varying betas
for (b in c(1:9)/10){
  #Set Beta
  Beta = b
  print(Beta)
  
  #Train
  rfmodel = train(pulsar~., 
                  machinefood, 
                  method = 'Rborist', 
                  trControl = trainControl(method = 'boot', summaryFunction = f1),
                  metric = 'F1')
  
  knnmodel = train(pulsar~., 
                   machinefood, 
                   method = 'knn', 
                   trControl = trainControl(method = 'boot', summaryFunction = f1), 
                   metric = 'F1')
  
  glmmodel = train(pulsar~., 
                   machinefood, 
                   method = 'glm', 
                   trControl = trainControl(method = 'boot', summaryFunction = f1),
                   metric = 'F1')
  
  nbmodel = train(pulsar~., 
                  machinefood, 
                  method = 'nb', 
                  trControl = trainControl(method = 'boot', summaryFunction = f1),
                  metric = 'F1')
  
  #Estimate performance
  temp = modelstats(rfmodel, Beta) %>% data.frame() %>% mutate(model = 'Rborist')
  allstats = union(allstats,temp)
  
  temp = modelstats(knnmodel, Beta) %>% data.frame() %>% mutate(model = 'knn')
  allstats = union(allstats, temp)
  
  temp = modelstats(glmmodel, Beta) %>% data.frame() %>% mutate(model = 'glm')
  allstats = union(allstats, temp)
  
  temp = modelstats(nbmodel, Beta) %>% data.frame() %>% mutate(model = 'nb')
  allstats = union(allstats, temp)
}

#Train models with specificty as metric
rfmodel = train(pulsar~., 
                machinefood, 
                method = 'Rborist', 
                trControl = trainControl(method = 'boot', summaryFunction = spec),
                metric = 'Specificity')

knnmodel = train(pulsar~., 
                machinefood, 
                method = 'knn', 
                trControl = trainControl(method = 'boot', summaryFunction = spec),
                metric = 'Specificity')

glmmodel = train(pulsar~., 
                machinefood, 
                method = 'glm', 
                trControl = trainControl(method = 'boot', summaryFunction = spec),
                metric = 'Specificity')

nbmodel = train(pulsar~., 
                machinefood, 
                method = 'nb', 
                trControl = trainControl(method = 'boot', summaryFunction = spec),
                metric = 'Specificity')
rfmodel = train(pulsar~., 
                machinefood, 
                method = 'Rborist', 
                trControl = trainControl(method = 'boot', summaryFunction = spec),
                metric = 'Specificity')

#Estimate performance
temp = modelstats(rfmodel, 1) %>% data.frame %>% mutate(model = 'Rborist Specificity')
allstats = union(allstats, temp)

temp = modelstats(knnmodel, 1) %>% data.frame %>% mutate(model = 'knn Specificity')
allstats = union(allstats, temp)

temp = modelstats(glmmodel, 1) %>% data.frame %>% mutate(model = 'glm Specificity')
allstats = union(allstats, temp)

temp = modelstats(nbmodel, 1) %>% data.frame %>% mutate(model = 'nb Specificity')
allstats = union(allstats, temp)


#Train models with sensitivity as metric
rfmodel = train(pulsar~., 
                machinefood, 
                method = 'Rborist', 
                trControl = trainControl(method = 'boot', summaryFunction = sens),
                metric = 'Sensitivity')

knnmodel = train(pulsar~., 
                 machinefood, 
                 method = 'knn', 
                 trControl = trainControl(method = 'boot', summaryFunction = sens),
                 metric = 'Sensitivity')

glmmodel = train(pulsar~., 
                 machinefood, 
                 method = 'glm', 
                 trControl = trainControl(method = 'boot', summaryFunction = sens),
                 metric = 'Sensitivity')

nbmodel = train(pulsar~., 
                machinefood, 
                method = 'nb', 
                trControl = trainControl(method = 'boot', summaryFunction = sens),
                metric = 'Sensitivity')
rfmodel = train(pulsar~., 
                machinefood, 
                method = 'Rborist', 
                trControl = trainControl(method = 'boot', summaryFunction = sens),
                metric = 'Sensitivity')

#Estimate performance
temp = modelstats(rfmodel, 1) %>% data.frame %>% mutate(model = 'Rborist Sensitivity')
allstats = union(allstats, temp)

temp = modelstats(knnmodel, 1) %>% data.frame %>% mutate(model = 'knn Sensitivity')
allstats = union(allstats, temp)

temp = modelstats(glmmodel, 1) %>% data.frame %>% mutate(model = 'glm Sensitivity')
allstats = union(allstats, temp)

temp = modelstats(nbmodel, 1) %>% data.frame %>% mutate(model = 'nb Sensitivity')
allstats = union(allstats, temp)




#Changes rownames (for aesthetics) and saves
rownames(allstats) = c(1:as.numeric(count(allstats)))
saveRDS(allstats, file = 'allstats.Rda')

#beta = .8 had best metrics so train with that and apply to test set
Beta = .8
rfmodel = train(pulsar~., 
                machinefood, 
                method = 'Rborist', 
                trControl = trainControl(method = 'boot', summaryFunction = f1),
                metric = 'F1')

temp = test %>% mutate(pulsarhat = predict(rfmodel, newdata = test))
ans = confusionMatrix(temp$pulsarhat, reference = temp$pulsar)
saveRDS(ans, 'ans.Rda')
# saveRDS(as.numeric(ans$byClass['Specificity']), 'Specificity.Rda')
# saveRDS(as.numeric(ans$overall['Accuracy']), 'Accuracy.Rda')

