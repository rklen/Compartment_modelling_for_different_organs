#functions
linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
fit_model<-function(param,k_number,input,j){
  k1<-abs(param[1])
  k2<-abs(param[2])
  if(k_number>=3){k3<-abs(param[3])}else{k3<-0}
  if(k_number==4){k4<-abs(param[4])}else{k4<-0}
  if(length(param)>=k_number+1){
    Va<-1/(1+exp(-param[k_number+1]))
    }else{Va<-0}
  if(length(param)==k_number+2){
    alpha<-1/(1+exp(-param[k_number+2]))
  }else{alpha<-1}
  v1<-0
  v2<-0
  for(i in 2:length(input)){
    if(i-j-1<1){c0<-0}else{
      if(i-j-1>length(input)){c0<-input[length(input)]}else{c0<-input[i-j-1]}
    }
    v1<-c(v1,k1*c0+(1-k2-k3)*v1[i-1]+k4*v2[i-1])
    v2<-c(v2,k3*v1[i-1]+(1-k4)*v2[i-1])
  }
  cT<-v1+v2
  v<-alpha*cT+Va*input
  return(v)
}
errorfunction<-function(param,k_number,input,organCurve,j){
  return(sum((organCurve-fit_model(param,k_number,input,j))^2))
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
optimizationAlg<-function(initialValues,k_number,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,
                k_number=k_number,input=input,
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
                '0046','0047','0048','0049','0050','0051','0052',
                '0053','0054','0055','0056','0057','0058','0059')
organNames<-c('Aorta','Brain','Myocardium','LLL','RLL','RML','LUL',
              'RUL','Liver','Spleen','Pancreas','Kidney, L','Kidney, R',
              'Colon','Bladder','Glu.max., R','Glu.med., L','Iliopsoas, R',
              'Humerus, R','10th rib, R','T5 vertebra')
modelNames<-c('1TCM','i2TCM','r2TCM',
              '1TCM+TDC','i2TCM+TDC','r2TCM+TDC',
              '1TCM+V_a','i2TCM+V_a','r2TCM+V_a',
              '1TCM+V_a+alpha','i2TCM+V_a+alpha','r2TCM+V_a+alpha',
              '1TCM+TDC+V_a','i2TCM+TDC+V_a','r2TCM+TDC+V_a',
              '1TCM+TDC+V_a+alpha','i2TCM+TDC+V_a+alpha','r2TCM+TDC+V_a+alpha')

#Models with no corrections
j<-0
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01)}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001)}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001)}
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:length(param)]<-param
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number-1,'.csv',sep=''),row.names=F)
  }
}

#Models with time delay correction only
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01)}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001)}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001)}
      errorsForj<-c()
      for(j in -30:30){
        param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
        err<-errorfunction(param,k_number,input,organCurve,j)
        errorsForj<-c(errorsForj,err)
      }
      j<-c(-30:30)[which.min(errorsForj)]
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:k_number]<-param
      df1[i,5]<-j
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number+2,'.csv',sep=''),row.names=F)
  }
}

#Models with V_a only
j<-0
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01,-log(9))}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001,-log(9))}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001,-log(9))}
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:k_number]<-param[1:k_number]
      df1[i,6]<-param[k_number+1]
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number+5,'.csv',sep=''),row.names=F)
  }
}

#Models with V_a and alpha only
j<-0
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01,-log(9),log(9))}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001,-log(9),log(9))}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001,-log(9),log(9))}
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:k_number]<-param[1:k_number]
      df1[i,6]<-param[k_number+1]
      df1[i,7]<-param[k_number+2]
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number+8,'.csv',sep=''),row.names=F)
  }
}

#Models with TDF and V_a
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01,-log(9))}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001,-log(9))}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001,-log(9))}
      errorsForj<-c()
      for(j in -30:30){
        param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
        err<-errorfunction(param,k_number,input,organCurve,j)
        errorsForj<-c(errorsForj,err)
      }
      j<-c(-30:30)[which.min(errorsForj)]
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:k_number]<-param[1:k_number]
      df1[i,5]<-j
      df1[i,6]<-param[k_number+1]
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number+11,'.csv',sep=''),row.names=F)
  }
}

#Models with TDF, V_a, and alpha
for(organIndex in 2:21){
  for(k_number in 2:4){
    df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=9)
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
      if(k_number==2){initialValues<-c(0.01,0.01,-log(9),log(9))}
      if(k_number==3){initialValues<-c(0.01,0.01,0.001,-log(9),log(9))}
      if(k_number==4){initialValues<-c(0.01,0.01,0.001,0.001,-log(9),log(9))}
      errorsForj<-c()
      for(j in -30:30){
        param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
        err<-errorfunction(param,k_number,input,organCurve,j)
        errorsForj<-c(errorsForj,err)
      }
      j<-c(-30:30)[which.min(errorsForj)]
      param<-optimizationAlg(initialValues,k_number,input,organCurve,j)
      points(times,fit_model(param,k_number,input,j),type='l')
      err<-errorfunction(param,k_number,input,organCurve,j)
      df1[i,1:k_number]<-param[1:k_number]
      df1[i,5]<-j
      df1[i,6]<-param[k_number+1]
      df1[i,7]<-param[k_number+2]
      df1[i,8]<-err/(length(times)-1)
      df1[i,9]<-meanRelError(fit_model(param,k_number,input,j),organCurve)
    }
    write.csv(df1,file=paste(organIndex,'_',k_number+14,'.csv',sep=''),row.names=F)
  }
}

bestModels<-c()
for(organIndex in 2:21){
  aics<-matrix(data=NA,nrow=length(studyNumbers)-1,ncol=18)
  for(l in 1:18){
    df1<-read.csv(paste(organIndex,'_',l,'.csv',sep=''))
    aics[,l]<-280*log(df1[c(1:8,10:56),8])+2*c(2,3,4,3,4,5,3,4,5,4,5,6,4,5,6,5,6,7)[l]
  }
  n<-which.min(colMedians(aics))
  bestModels<-c(bestModels,n)
}

string<-''
for(organIndex in 2:21){
  n<-bestModels[organIndex-1]
  string<-c(string,organNames[organIndex],'&',modelNames[n],'&')
  if(n==5){
    string<-c(string,wilcoxSymbol(aics[,5],aics[,4]),'&')
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,5],aics[,6]),'&')
    string<-c(string,wilcoxSymbol(aics[,5],aics[,2]),'&')
    string<-c(string,'-&-newRow')
  }
  if(n==6){
    string<-c(string,wilcoxSymbol(aics[,6],aics[,4]),'&')
    string<-c(string,wilcoxSymbol(aics[,6],aics[,5]),'&')
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,6],aics[,3]),'&')
    string<-c(string,'-&-newRow')
  }
  if(n==13){
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,13],aics[,14]),'&')
    string<-c(string,wilcoxSymbol(aics[,13],aics[,15]),'&')
    string<-c(string,wilcoxSymbol(aics[,13],aics[,7]),'&')
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,13],aics[,4]),'newRow')
  }
  if(n==16){
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,16],aics[,17]),'&')
    string<-c(string,wilcoxSymbol(aics[,16],aics[,18]),'&')
    string<-c(string,wilcoxSymbol(aics[,16],aics[,10]),'&')
    string<-c(string,wilcoxSymbol(aics[,16],aics[,13]),'&')
    string<-c(string,wilcoxSymbol(aics[,16],aics[,4]),'newRow')
  }
  if(n==18){
    string<-c(string,wilcoxSymbol(aics[,18],aics[,16]),'&')
    string<-c(string,wilcoxSymbol(aics[,18],aics[,17]),'&')
    string<-c(string,'-&')
    string<-c(string,wilcoxSymbol(aics[,18],aics[,12]),'&')
    string<-c(string,wilcoxSymbol(aics[,18],aics[,15]),'&')
    string<-c(string,wilcoxSymbol(aics[,13],aics[,6]),'newRow')
  }
  organIndex<-organIndex+1
}
cat(string)

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  df1<-read.csv(paste(organIndex,'_',bestModels[organIndex-1],'.csv',sep=''))
  for(l in 1:7){
    if(l<=4){
      string<-c(string,'&',signif(mean(abs(df1[c(1:8,10:56),l]))*60,3),'plusMinus',signif(sd(abs(df1[c(1:8,10:56),l]))*60,3))
    }
    if(l==5){
      string<-c(string,'&',signif(mean(df1[c(1:8,10:56),l]),3),'plusMinus',signif(sd(df1[c(1:8,10:56),l]),3))
    }
    if(l>=6){
      if(bestModels[organIndex-1]<7){
        string<-c(string,'&-')
      }else{
        string<-c(string,'&',signif(mean(1/(1+exp(-df1[c(1:8,10:56),l])))*100,3),'plusMinus',signif(sd(1/(1+exp(-df1[c(1:8,10:56),l])))*100,3))
      }
    }
  }
  string<-c(string,'newRow')
}
cat(string)

string<-''
for(organIndex in 2:21){
  string<-c(string,organNames[organIndex])
  df1<-read.csv(paste(organIndex,'_',bestModels[organIndex-1],'.csv',sep=''))
  mres1<-df1[c(1:8,10:56),9]
  for(index in 2:21){
    if(index!=organIndex){
      df<-read.csv(paste(index,'_',bestModels[index-1],'.csv',sep=''))
      mres<-df[c(1:8,10:56),9]
      string<-c(string,'&',wilcoxSymbol(mres1,mres))
    }else{
      string<-c(string,'&-')
    }
  }
  string<-c(string,'newRow')
}
cat(string)

organIndex<-21
number<-'0001'
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
if(bestModels[organIndex-1]==5){
  param<-c(df1[1,1],df1[1,2],df1[1,3])
  j<-df1[1,5]
  k_number<-3
}
if(bestModels[organIndex-1]==6){
  param<-c(df1[1,1],df1[1,2],df1[1,3],df1[1,4])
  j<-df1[1,5]
  k_number<-4
}
if(bestModels[organIndex-1]==10){
  param<-c(df1[1,1],df1[1,2],df1[1,6],df1[1,7])
  j<-df1[1,5]
  k_number<-2
}
if(bestModels[organIndex-1]==12){
  param<-c(df1[1,1],df1[1,2],df1[1,3],df1[1,4],df1[1,6],df1[1,7])
  j<-df1[1,5]
  k_number<-4
}
points(times,fit_model(param,k_number,input,j),type='l',col='blue',lwd=2)

