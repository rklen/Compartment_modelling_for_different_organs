linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
fit_1Tmodel<-function(k_val,input){
  v1<-0
  for(k in 2:length(input)){
    v1<-c(v1,v1[k-1]+input[k-1]*abs(k_val[1])-v1[k-1]*abs(k_val[2]))
  }
  return(v1)
}
fit_2Tmodel<-function(k_val,input){
  v1<-0
  v2<-0
  for(k in 2:length(input)){
    v1<-c(v1,v1[k-1]+input[k-1]*abs(k_val[1])-v1[k-1]*(abs(k_val[2])+abs(k_val[3])))
    v2<-c(v2,v2[k-1]+v1[k-1]*abs(k_val[3]))
  }
  v<-v1+v2
  return(v)
}
errorfunction<-function(k_val,T_number,input,organCurve){
  if(T_number==1){
    return(sum((organCurve-fit_1Tmodel(k_val,input))^2))
  }
  if(T_number==2){
    return(sum((organCurve-fit_2Tmodel(k_val,input))^2))
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
optimizationAlg<-function(initialValues,T_number,input,organCurve){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,T_number=T_number,input=input,
           organCurve=organCurve,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
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
    k<-optimizationAlg(initialValues,T_number,input,organCurve)
    k<-c(abs(k[1]),abs(k[2]))
    df1[i,3]<-k[1]
    df1[i,4]<-k[2]
    points(times,fit_1Tmodel(k,input),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve)
    df1[i,5]<-err/(length(times)-1)
    df1[i,6]<-meanRelError(fit_1Tmodel(k,input),organCurve)
    T_number<-2
    initialValues<-c(0.01,0.01,0.001)
    k<-optimizationAlg(initialValues,T_number,input,organCurve)
    k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
    df1[i,7]<-k[1]
    df1[i,8]<-k[2]
    df1[i,9]<-k[3]
    points(times,fit_2Tmodel(k,input),type='l')
    err<-errorfunction(k,T_number=T_number,input,organCurve)
    df1[i,10]<-err/(length(times)-1)
    df1[i,11]<-meanRelError(fit_2Tmodel(k,input),organCurve)
  }
  write.csv(df1,file=paste('lok',organIndex,'.csv',sep=''),row.names=F)
}

organIndex<-21
df1<-read.csv(paste('lok',organIndex,'.csv',sep=''))
v<-c()
for(i in c(1,3:11)){
  v1<-median(df1[,i])
  if(v1>1){
    v1<-round(v1,2)
  }
  if(v1<1 & v1>0.1){
    v1<-round(v1,3)
  }
  if(v1<0.1 & v1>0.01){
    v1<-round(v1,4)
  }
  if(v1<0.01 & v1>0.001){
    v1<-round(v1,5)
  }
  v<-c(v,v1)
} 
#print(v)
aicsFor1T<-280*log(df1[,5])+2*2
aicsFor2T<-280*log(df1[,10])+2*3
(median(df1[,10])<median(df1[,5]))
wilcox.test(aicsFor1T,aicsFor2T,paired=TRUE)$p.value
(median(df1[,11])<median(df1[,6]))
wilcox.test(aicsFor1T,aicsFor2T,paired=TRUE)$p.value
(median(aicsFor2T)<median(aicsFor1T))
wilcox.test(aicsFor1T,aicsFor2T,paired=TRUE)$p.value

organNames<-c('aorta','brain','myocard','lll','rll','rml','lul',
              'rul','liver','spleen','pancreas','lkidney','rkidney',
              'colon','bladder','rglumax','lglumed','rilio',
              'rhum','10rib','T5')
modelIndex<-c(NA,2,2,2,2,2,1,1,1,1,1,2,2,2,1,1,2,2,2,2,2)
modelIndex<-rep(2,21)
organIndex<-21
df<-read.csv(paste('lok',organIndex,'.csv',sep=''))
if(modelIndex[organIndex]==1){
  errors<-df[,6]
}else{
  errors<-df[,11]
}
results<-c()
for(index in 2:21){
  if(index==organIndex){
    s<-'& -'
  }else{
    df1<-read.csv(paste('lok',index,'.csv',sep=''))
    if(modelIndex[index]==1){
      errors1<-df1[,6]
    }else{
      errors1<-df1[,11]
    }
    s<-'&'
    if(median(errors)<median(errors1)){
      if(wilcox.test(errors,errors1,paired=TRUE)$p.value<0.05){
        s<-'& *'
      }
      if(wilcox.test(errors,errors1,paired=TRUE)$p.value<0.01){
        s<-'& **'
      }
      if(wilcox.test(errors,errors1,paired=TRUE)$p.value<0.001){
        s<-'& ***'
      }
    }
  }
  results<-c(results,s)
}
print(results,quote=FALSE)

df<-read.csv('C:/Users/oonar/Downloads/array_koveri0001.csv')
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,21)
input<-c()
for(j in 1:length(times)){
  input[j]<-linearInter(times[j],t,df[1,])
}
plot(times,input,type='l',xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5)

organIndex<-7
organCurve<-c()
for(j in 1:length(times)){
  organCurve[j]<-linearInter(times[j],t,df[organIndex,])
}
plot(times,organCurve,type='l',xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5,
     #ylim=c(0,4.7)
     )
T_number<-1
initialValues<-c(0.01,0.01)
k<-optimizationAlg(initialValues,T_number,input,organCurve)
k<-c(abs(k[1]),abs(k[2]))
print(k)
points(times,fit_1Tmodel(k,input),type='l',col='gray',lwd=2)
err<-errorfunction(k,T_number=T_number,input,organCurve)
print(err/(length(times)-1))
print(meanRelError(fit_1Tmodel(k,input),organCurve))
T_number<-2
initialValues<-c(0.01,0.01,0.001)
k<-optimizationAlg(initialValues,T_number,input,organCurve)
k<-c(abs(k[1]),abs(k[2]),abs(k[3]))
print(k)
points(times,fit_2Tmodel(k,input),type='l',col='blue',lwd=2)
err<-errorfunction(k,T_number=T_number,input,organCurve)
print(err/(length(times)-1))
print(meanRelError(fit_2Tmodel(k,input),organCurve))

x<-c()
y<-c()
z<-c()
for(organIndex in 2:21){
  df1<-read.csv(paste('lok',organIndex,'.csv',sep=''))
  x<-c(x,median(df1[,2]))
  y<-c(y,median(df1[,6]))
  z<-c(z,median(df1[,11]))
}
plot(x,y,xlab='',ylab='',lwd=2,cex.axis=1.5,cex.lab=1.5,
     ylim=c(0,0.25))
i=1
organNames[i+1]
points(x[i],y[i],pch=3)
cor(x,y)
cor(x,y,method='spearman')
plot(x,z)
cor(x,z)
cor(x,z,method='spearman')