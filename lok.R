#functions
linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
fit_model<-function(param,input,j){
  k1<-abs(param[1])
  k2<-abs(param[2])
  if(length(param)>=3){Va<-1/(1+exp(-param[3]))}else{Va<-0}
  if(length(param)==4){alpha<-1/(1+exp(-param[4]))}else{alpha<-1}
  cT<-0
  for(i in 2:length(input)){
    if(i-j-1<1){c0<-0}else{
      if(i-j-1>length(input)){c0<-input[length(input)]}else{c0<-input[i-j-1]}
    }
    cT<-c(cT,k1*c0+(1-k2)*cT[i-1])
  }
  v<-alpha*(1-Va)*cT+Va*input
  return(v)
}
errorfunction<-function(param,input,organCurve,j){
  return(sum((organCurve-fit_model(param,input,j))^2))
}
meanRelError<-function(modelCurve,organCurve){
  r<-0
  for(i in 2:length(organCurve)){
    if(modelCurve[i]!=organCurve[i]){
      r<-r+abs(modelCurve[i]-organCurve[i])/organCurve[i]
    }
  }
  r<-r/(length(organCurve)-1)
  return(r)
}
optimizationAlg<-function(initialValues,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
           organCurve=organCurve,j=j,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
colMedians<-function(d){
  medians<-c()
  for(i in 1:ncol(d)){
    medians[i]<-median(d[,i])
  }
  return(medians)
}
wilcoxSymbol<-function(x,y){
  s<-''
  p<-wilcox.test(x,y,alternative='less',paired=TRUE)$p.value
  if(p<0.05){s<-'*'}
  if(p<0.01){s<-'**'}
  if(p<0.001){s<-'***'}
  return(s)
}

t<-c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,
     160,190,220,250,280)
times<-c(0:280)
studyNumbers<-c('0001','0002','0003','0004','0005','0006','0007',
                '0008','0009','0010','0011','0012','0013','0014',
                '0015','0016','0017','0018','0019','0020','0021',
                '0022','0023','0024','0025','0026','0027','0029',
                '0030','0031','0032','0033','0035','0036','0037',
                '0038','0039','0040','0041','0042','0043','0044',
                '0045','0046','0047','0048','0049','0050','0051',
                '0052','0053','0054','0055','0056','0057','0058',
                '0059','0060')
organNames<-c('Aorta','Brain','Myocardium','LLL','RLL','RML','LUL',
              'RUL','Liver','Spleen','Pancreas','Kidney, L','Kidney, R',
              'Colon','Bladder','Glu.max., R','Glu.med., L','Iliopsoas, R',
              'Humerus, R','10th rib, R','T5 vertebra')
modelNames<-c('1TCM','1TCM+TDC','1TCM+V_a','1TCM+V_a+alpha',
              '1TCM+TDC+V_a','1TCM+TDC+V_a+alpha')

#1TCM with no corrections
j<-0
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01)
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',1,'.csv',sep=''),row.names=F)
}

#1TCM with time delay correction only
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01)
    errorsForj<-c()
    for(j in -30:30){
      param<-optimizationAlg(initialValues,input,organCurve,j)
      err<-errorfunction(param,input,organCurve,j)
      errorsForj<-c(errorsForj,err)
    }
    j<-c(-30:30)[which.min(errorsForj)]
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,5]<-j
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',2,'.csv',sep=''),row.names=F)
}

#1TCM with V_a only
j<-0
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01,-log(9))
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',3,'.csv',sep=''),row.names=F)
}

#1TCM with V_a and alpha only
j<-0
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01,-log(9),log(9))
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',4,'.csv',sep=''),row.names=F)
}

#1TCM with TDF and V_a
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01,-log(9))
    errorsForj<-c()
    for(j in -30:30){
      param<-optimizationAlg(initialValues,input,organCurve,j)
      err<-errorfunction(param,input,organCurve,j)
      errorsForj<-c(errorsForj,err)
    }
    j<-c(-30:30)[which.min(errorsForj)]
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,5]<-j
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',5,'.csv',sep=''),row.names=F)
}

#1TCM with TDF, V_a, and alpha
for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(l in 1:length(times)){
      input[l]<-linearInter(times[l],t,df[1,])
      organCurve[l]<-linearInter(times[l],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    initialValues<-c(0.01,0.01,-log(9),log(9))
    errorsForj<-c()
    for(j in -30:30){
      param<-optimizationAlg(initialValues,input,organCurve,j)
      err<-errorfunction(param,input,organCurve,j)
      errorsForj<-c(errorsForj,err)
    }
    j<-c(-30:30)[which.min(errorsForj)]
    param<-optimizationAlg(initialValues,input,organCurve,j)
    points(times,fit_model(param,input,j),type='l')
    err<-errorfunction(param,input,organCurve,j)
    df1[i,1:length(param)]<-param
    df1[i,5]<-j
    df1[i,6]<-err/(length(times)-1)
    df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
  }
  write.csv(df1,file=paste(organIndex,'_',6,'.csv',sep=''),row.names=F)
}

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  df1<-read.csv(paste(organIndex,'_',6,'.csv',sep=''))
  string<-c(string,'&',signif(mean(abs(df1[,1]))*60,3),'plusMinus',signif(sd(abs(df1[,1]))*60,3))
  string<-c(string,'&',signif(mean(abs(df1[,2]))*60,3),'plusMinus',signif(sd(abs(df1[,2]))*60,3))
  string<-c(string,'&',signif(mean(abs(df1[,1]/df1[,2])),3),'plusMinus',signif(sd(abs(df1[,1]/df1[,2])),3))
  string<-c(string,'&',signif(mean(1/(1+exp(-df1[,3])))*100,3),'plusMinus',signif(sd(1/(1+exp(-df1[,3])))*100,3))
  string<-c(string,'&',signif(mean(1/(1+exp(-df1[,4])))*100,3),'plusMinus',signif(sd(1/(1+exp(-df1[,4])))*100,3))
  string<-c(string,'&',signif(mean(df1[,5]),3),'plusMinus',signif(sd(df1[,5]),3))
  string<-c(string,'newRow')
}
cat(string)

bestModels<-c()
for(organIndex in 2:21){
  aics<-matrix(data=NA,nrow=length(studyNumbers),ncol=6)
  for(l in 1:6){
    df1<-read.csv(paste(organIndex,'_',l,'.csv',sep=''))
    aics[,l]<-280*log(df1[,6])+2*c(2,3,3,4,4,5)[l]
  }
  n<-which.min(colMedians(aics))
  bestModels<-c(bestModels,n)
}

string<-''
for(organIndex in 2:21){
  n<-bestModels[organIndex-1]
  aics<-matrix(data=NA,nrow=length(studyNumbers),ncol=6)
  for(l in 1:6){
    df1<-read.csv(paste(organIndex,'_',l,'.csv',sep=''))
    aics[,l]<-280*log(df1[,6])+2*c(2,3,3,4,4,5)[l]
  }
  string<-c(string,organNames[organIndex-1],'&',modelNames[n])
  for(l in c(1,3,4,2,5,6)){
    if(l!=n){
      s<-wilcoxSymbol(aics[,n],aics[,l])
    }else{s<-'-'}
    string<-c(string,'&',s)
  }
  string<-c(string,'newRow ')
}
cat(string)

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  n<-bestModels[organIndex-1]
  df1<-read.csv(paste(organIndex,'_',n,'.csv',sep=''))
  string<-c(string,'&',signif(median(abs(df1[,1]))*60,3))
  string<-c(string,'&',signif(median(abs(df1[,2]))*60,3))
  string<-c(string,'&',signif(median(abs(df1[,1]/df1[,2])),3))
  string<-c(string,'&',signif(median(1/(1+exp(-df1[,3])))*100,3))
  string<-c(string,'&',signif(median(1/(1+exp(-df1[,4])))*100,3))
  string<-c(string,'&',signif(median(df1[,5]),3))
  string<-c(string,'newRow')
}
cat(string)

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  n<-bestModels[organIndex-1]
  df1<-read.csv(paste(organIndex,'_',n,'.csv',sep=''))
  string<-c(string,'&',signif(median(df1[,6]),3))
  string<-c(string,'&',signif(median(df1[,7]),3))
  string<-c(string,'&',signif(median(280*log(df1[,6])+2*c(2,3,3,4,4,5)[n],3)))
  string<-c(string,'newRow')
}
cat(string)

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  df1<-read.csv(paste(organIndex,'_',bestModels[organIndex-1],'.csv',sep=''))
  mres1<-df1[,7]
  for(index in 2:21){
    if(index!=organIndex){
      df<-read.csv(paste(index,'_',bestModels[index-1],'.csv',sep=''))
      mres<-df[,7]
      string<-c(string,'&',wilcoxSymbol(mres1,mres))
    }else{
      string<-c(string,'&-')
    }
  }
  string<-c(string,'newRow')
}
cat(string)

df<-read.table('C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt')
df<-df[2:60,]
w<-c(90,94,79,85,73,82,81,108,63,92,64,117,99,59,89,100,103,101,78,86,87,77,114,81,90,78,65,94,121,61,95,64,60,92,101,69,62,64,132,
     98,130,67,123,80,NA,102,93,76,91,72,87,88,78,113,110,81,95,69,107,NA)
dose<-c(355,356,351,308,316,336,364,354,350,406,357,319,334,350,356,341,339,332,324,322,408,402,364,385,363,334,380,364,336,351,
        362,366,396,390,359,364,334,345,379,390,348,366,321,332,295,370,356,371,323,307,374,383,333,381,349,338,353,360,362,366)

organIndex<-1
number<-'0014'
filepath<-paste('C:/Users/Oona/Documents/Tpc/lok/array_koveri',number,'.csv',sep='')
df<-read.csv(filepath)
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,21)
input<-c()
organCurve<-c()
for(l in 1:length(times)){
  input[l]<-linearInter(times[l],t,df[1,])
  organCurve[l]<-linearInter(times[l],t,df[organIndex,])
}
plot(times,organCurve,type='l',xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5,
     #ylim=c(0,7.5)
)
df1<-read.csv(paste(organIndex,'_',bestModels[organIndex-1],'.csv',sep=''))
if(is.nan(df1[1,3])){
  param<-c(df1[1,1],df1[1,2])
}else{
  if(is.nan(df1[1,4])){
    param<-c(df1[1,1],df1[1,2],df1[1,3])
  }else{
    param<-c(df1[1,1],df1[1,2],df1[1,3],df1[1,4])
  }
}
j<-df1[1,5]
points(times,fit_model(param,input,j),type='l',lwd=2,col='blue')

filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df<-read.table(filepath)
df<-df[c(2:28,30:34,36:61),]
sum(df$V3=='Female')
sum(df$V3=='Male')
length(studyNumbers)
mean(as.numeric(df$V2))
sd(as.numeric(df$V2))
min(as.numeric(df$V2))
max(as.numeric(df$V2))
