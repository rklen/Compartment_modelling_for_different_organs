linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
fit_1Tmodel<-function(k_val,input,timeDelay){
  v1<-rep(0,timeDelay+1)
  for(k in (timeDelay+2):length(input)){
    v1<-c(v1,v1[k-1]+input[k-1-timeDelay]*abs(k_val[1])-v1[k-1]*abs(k_val[2]))
  }
  return(v1)
}
fit_2Tmodel<-function(k_val,input,timeDelay){
  v1<-rep(0,timeDelay+1)
  v2<-rep(0,timeDelay+1)
  for(k in (timeDelay+2):length(input)){
    v1<-c(v1,v1[k-1]+input[k-1-timeDelay]*abs(k_val[1])-v1[k-1]*(abs(k_val[2])+abs(k_val[3])))
    v2<-c(v2,v2[k-1]+v1[k-1]*abs(k_val[3]))
  }
  v<-v1+v2
  return(v)
}
errorfunction<-function(k_val,T_number,input,organCurve,timeDelay){
  if(T_number==1){
    return(sum((organCurve-fit_1Tmodel(k_val,input,timeDelay))^2))
  }
  if(T_number==2){
    return(sum((organCurve-fit_2Tmodel(k_val,input,timeDelay))^2))
  }
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
optimizationAlg<-function(initialValues,T_number,input,organCurve,timeDelay){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,T_number=T_number,input=input,
           organCurve=organCurve,timeDelay=timeDelay,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
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
                '0046','0047')
organNames<-c('aorta','brain','myocard','lll','rll','rml','lul',
              'rul','liver','spleen','pancreas','lkidney','rkidney',
              'colon','bladder','rglumax','lglumed','rilio',
              'rhum','10rib','T5')

for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=22)
  for(timeDelay in 0:10){
    for(i in 1:length(studyNumbers)){
      number<-studyNumbers[i]
      filepath<-paste('C:/Users/oonar/Downloads/array_koveri',
                      number,'.csv',sep='')
      df<-read.csv(filepath)
      df<-as.matrix(df)
      colnames(df)<-NULL
      rownames(df)<-NULL
      df[,1]<-rep(0,21)
      input<-c()
      organCurve<-c()
      for(j in 1:length(times)){
        input[j]<-linearInter(times[j],t,df[1,])
        organCurve[j]<-linearInter(times[j],t,df[organIndex,])
      }
      #plot(times,input,type='l')
      plot(times,organCurve,type='l')
      T_number<-1
      initialValues<-c(0.01,0.01)
      k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
      k<-c(abs(k[1]),abs(k[2]))
      points(times,fit_1Tmodel(k,input,timeDelay),type='l')
      err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
      df1[i,(timeDelay+1)]<-err/(length(times)-1)
      T_number<-2
      initialValues<-c(0.01,0.01,0.001)
      k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
      k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
      points(times,fit_2Tmodel(k,input,timeDelay),type='l')
      err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
      df1[i,(timeDelay+12)]<-err/(length(times)-1)
    }
  }
  write.csv(df1,file=paste('tdlok',organIndex,'.csv',sep=''),row.names=F)
}

tdc_df<-matrix(data=NA,nrow=length(studyNumbers),ncol=40)
for(organIndex in 2:21){
  df1<-read.csv(paste('tdlok',organIndex,'.csv',sep=''))
  v0<-c()
  v1<-c()
  for(i in 1:44){
    tdc_df[i,organIndex-1]<-which.min(df1[i,1:11])-1
    tdc_df[i,organIndex+19]<-which.min(df1[i,12:22])-1
  }
}
write.csv(tdc_df,file='tdc_df.csv',row.names=F)

for(organIndex in 2:21){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=11)
  colnames(df1)<-c('maxSUV','loc','1T_k1','1T_k2','1T_mse','1T_mre',
                   '2T_k1','2T_k2','2T_k3','2T_mse','2T_mre')
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/oonar/Downloads/array_koveri',
                    number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(j in 1:length(times)){
      input[j]<-linearInter(times[j],t,df[1,])
      organCurve[j]<-linearInter(times[j],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    df1[i,1]<-max(organCurve)
    df1[i,2]<-which.max(organCurve)-1
    T_number<-1
    initialValues<-c(0.01,0.01)
    timeDelay<-tdc_df[i,organIndex-1]
    k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
    k<-c(abs(k[1]),abs(k[2]))
    df1[i,3]<-k[1]
    df1[i,4]<-k[2]
    points(times,fit_1Tmodel(k,input,timeDelay),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
    df1[i,5]<-err/(length(times)-1)
    df1[i,6]<-meanRelError(fit_1Tmodel(k,input,timeDelay),organCurve)
    T_number<-2
    initialValues<-c(0.01,0.01,0.001)
    timeDelay<-tdc_df[i,organIndex+19]
    k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
    k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
    df1[i,7]<-k[1]
    df1[i,8]<-k[2]
    df1[i,9]<-k[3]
    points(times,fit_2Tmodel(k,input,timeDelay),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
    df1[i,10]<-err/(length(times)-1)
    df1[i,11]<-meanRelError(fit_2Tmodel(k,input,timeDelay),organCurve)
  }
  write.csv(df1,file=paste('lok',organIndex,'.csv',sep=''),row.names=F)
}

organIndex<-21
df1<-read.csv(paste('lok',organIndex,'.csv',sep=''))
string<-''
for(i in c(1,2,3,4,7,8,9)){
  string<-paste(string,' & ',signif(mean(df1[,i]),3),sep='')
  if(i==4){
    string<-paste(string,' & ',signif(mean(tdc_df[,organIndex-1]),3),sep='')
  }
  if(i==9){
    string<-paste(string,' & ',signif(mean(tdc_df[,organIndex+19]),3),sep='')
  }
}
print(string)
string<-''
for(i in c(1,2,3,4,7,8,9)){
  string<-paste(string,' & $pm$',signif(sd(df1[,i]),3),sep='')
  if(i==4){
    string<-paste(string,' & $pm$',signif(sd(tdc_df[,organIndex-1]),3),sep='')
  }
  if(i==9){
    string<-paste(string,' & $pm$',signif(sd(tdc_df[,organIndex+19]),3),sep='')
  }
}
print(string)

organIndex<-21
df<-read.csv(paste('lok',organIndex,'.csv',sep=''))
aicsFor1T<-280*log(df[,5])+2*3
aicsFor2T<-280*log(df[,10])+2*4
df1<-read.csv(paste('tdlok',organIndex,'.csv',sep=''))
aicsFor1T_noTDC<-280*log(df1[,1])+2*2
aicsFor2T_noTDC<-280*log(df1[,12])+2*3
finModIndex<-which.min(c(
  median(aicsFor1T),median(aicsFor2T),median(aicsFor1T_noTDC),
  median(aicsFor2T_noTDC)))
string<-c('TDC T1','TDC T2','non-TDC T1','non-TDC T2')[finModIndex]
if(finModIndex==1){
  write.csv(df[,6],file=paste('mre',organIndex,'.csv',sep=''),row.names=F)
  string<-paste(string,' & ',signif(median(df[,5]),3),' & ',
                signif(median(df[,6]),3),sep='')
}
if(finModIndex==2){
  write.csv(df[,11],file=paste('mre',organIndex,'.csv',sep=''),row.names=F)
  string<-paste(string,' & ',signif(median(df[,10]),3),' & ',
                signif(median(df[,11]),3),sep='')
}
if(finModIndex==3){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=11)
  colnames(df1)<-c('maxSUV','loc','1T_k1','1T_k2','1T_mse','1T_mre',
                   '2T_k1','2T_k2','2T_k3','2T_mse','2T_mre')
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/oonar/Downloads/array_koveri',
                    number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(j in 1:length(times)){
      input[j]<-linearInter(times[j],t,df[1,])
      organCurve[j]<-linearInter(times[j],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    T_number<-1
    initialValues<-c(0.01,0.01)
    timeDelay<-0
    k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
    k<-c(abs(k[1]),abs(k[2]))
    points(times,fit_1Tmodel(k,input,timeDelay),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
    df1[i,5]<-err/(length(times)-1)
    df1[i,6]<-meanRelError(fit_1Tmodel(k,input,timeDelay),organCurve)
  }
  write.csv(df1[,6],file=paste('mre',organIndex,'.csv',sep=''),row.names=F)
  string<-paste(string,' & ',signif(median(df1[,5]),3),' & ',
                signif(median(df1[,6]),3),sep='')
}
if(finModIndex==4){
  df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=11)
  colnames(df1)<-c('maxSUV','loc','1T_k1','1T_k2','1T_mse','1T_mre',
                   '2T_k1','2T_k2','2T_k3','2T_mse','2T_mre')
  for(i in 1:length(studyNumbers)){
    number<-studyNumbers[i]
    filepath<-paste('C:/Users/oonar/Downloads/array_koveri',
                    number,'.csv',sep='')
    df<-read.csv(filepath)
    df<-as.matrix(df)
    colnames(df)<-NULL
    rownames(df)<-NULL
    df[,1]<-rep(0,21)
    input<-c()
    organCurve<-c()
    for(j in 1:length(times)){
      input[j]<-linearInter(times[j],t,df[1,])
      organCurve[j]<-linearInter(times[j],t,df[organIndex,])
    }
    #plot(times,input,type='l')
    plot(times,organCurve,type='l')
    T_number<-2
    initialValues<-c(0.01,0.01,0.001)
    timeDelay<-tdc_df[i,organIndex+19]
    k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
    k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
    points(times,fit_2Tmodel(k,input,timeDelay),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve,timeDelay)
    df1[i,10]<-err/(length(times)-1)
    df1[i,11]<-meanRelError(fit_2Tmodel(k,input,timeDelay),organCurve)
  }
  write.csv(df1[,11],file=paste('mre',organIndex,'.csv',sep=''),row.names=F)
  string<-paste(string,' & ',signif(median(df1[,10]),3),' & ',
                signif(median(df1[,11]),3),sep='')
}
string<-paste(string,' & ',wilcoxSymbol(aicsFor2T,aicsFor1T),sep='')
string<-paste(string,' & ',wilcoxSymbol(aicsFor2T_noTDC,aicsFor1T_noTDC),sep='')
string<-paste(string,' & ',wilcoxSymbol(aicsFor2T,aicsFor2T_noTDC),sep='')
string<-paste(string,' & ',wilcoxSymbol(aicsFor1T,aicsFor1T_noTDC),sep='')
print(string)

organIndex<-2
errors<-read.csv(paste('mre',organIndex,'.csv',sep=''))[,1]
string<-''
for(index in 2:21){
  if(index==organIndex){
    string<-paste(string,'& -',sep='')
  }else{
    errors1<-read.csv(paste('mre',index,'.csv',sep=''))[,1]
    string<-paste(string,' & ',wilcoxSymbol(errors,errors1),sep='')
  }
}
print(string)

i<-1
df<-read.csv('C:/Users/oonar/Downloads/array_koveri0001.csv')
df<-read.csv(filepath)
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,21)
input<-c()
for(j in 1:length(times)){
  input[j]<-linearInter(times[j],t,df[1,])
}
plot(times,input,type='l',xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5)
organIndex<-21
organCurve<-c()
for(j in 1:length(times)){
  organCurve[j]<-linearInter(times[j],t,df[organIndex,])
}
plot(times,organCurve,type='l',xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5,
     #ylim=c(0,7.5)
     )
T_number<-1
initialValues<-c(0.01,0.01)
timeDelay<-tdc_df[i,organIndex-1]
k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
k<-c(abs(k[1]),abs(k[2]))
points(times,fit_1Tmodel(k,input,timeDelay),type='l',col='gray',lwd=2)
T_number<-2
initialValues<-c(0.01,0.01,0.001)
timeDelay<-tdc_df[i,organIndex+19]
k<-optimizationAlg(initialValues,T_number,input,organCurve,timeDelay)
k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
points(times,fit_2Tmodel(k,input,timeDelay),type='l',col='blue',lwd=2)

x<-c()
y<-c()
for(organIndex in 2:21){
  df1<-read.csv(paste('lok',organIndex,'.csv',sep=''))
  x<-c(x,median(df1[,1]))
  y<-c(y,median(read.csv(paste('mre',organIndex,'.csv',sep=''))[,1]))
}
plot(x,y,xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5,
     ylim=c(0,0.20))
cor(x,y)
cor(x,y,method='spearman')
i=1
organNames[i+1]
points(x[i],y[i],pch=3)
