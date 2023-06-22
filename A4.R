matrix.2ndorder.make<-function(x, only.quad=F){
  x0<-x
  dimn<-dimnames(x)[[2]] #extract the names of the variables
  num.col<-length(x[1,]) # how many columns
  for(i in 1:num.col){
    # if we are doing all 2nd order
    if(!only.quad){
      for(j in i:num.col){
        x0<-cbind(x0,as.numeric(x[,i])*as.numeric(x[,j]))
        dimn<-c(dimn,paste(dimn[i],dimn[j],sep=""))
        #create interaction dimnames
      }
    }
    else{
      #in here only if doing only squared terms
      x0<-cbind(x0,as.numeric(x[,i])*as.numeric(x[,i]))
      dimn<-c(dimn,paste(dimn[i],"2",sep="")) # squared dimmension names
    }
  }
  dimnames(x0)[[2]]<-dimn
  x0
}
#sum of the absolute values of x
sumabs<-function(x){sum(abs(x))}
betanorm.lars<-function(str){
  v1<-apply(str$beta,1,sumabs)
  v1/max(v1)
}
regpluspress<-function(x,y){
  #least squares fit of x and y
  str<-lsfit(x,y)
  #calculates press statistic
  press<-sum((str$resid/(1-hat(x)))^2)
  #assigns press and coeff to attributes of str
  str$press<-press
  beta<-str$coefficients
  str$beta<-beta
  str
}
leaps.then.press<-function(xmat,yvec,ncheck=10,print.ls=F)
{
  #performs leaps on xmat and yvec 
  leaps.str<-leaps(xmat,yvec)
  # plots the size against corrosponding cps
  plot(leaps.str$size, log(leaps.str$Cp))
  #saves cpstatistics
  z1<-leaps.str$Cp
  #orders cp and saves the order
  o1<-order(z1)
  #saves ncheck models with respect to cp in order
  matwhich<-(leaps.str$which[o1,])[1:ncheck,]
  z2<-z1[o1][1:ncheck]
  #initializes press vec to use in following code
  pressvec<-NULL
  #calculates the press for the ncheck best models with regard to cp
  for(i in 1:ncheck){
    ls.str0<-regpluspress(xmat[,matwhich[i,]],yvec)
    print(i)
    print(paste("Press=",ls.str0$press))
    #saves each press value into press vec
    pressvec<-c(pressvec,ls.str0$press)
    parvec<-matwhich[i,]
    #number of parameters
    npar<-sum(parvec)
    #calculates the mean predicted square error
    print(paste("MPSE=",ls.str0$press/(length(yvec)-(npar+1))))
    print(paste("Cp=",z2[i]))
    if(print.ls){
      print(" Model Statistics ")
      ls.print(ls.str0)
    }
  }
  #saves the index of the model with the smallest press stat
  modin<-which.min(pressvec)
  ls.strB<-regpluspress(xmat[,matwhich[modin,]],yvec)
  print("===============")
  print("Best Model")
  print(paste("Press=",ls.strB$press))
  parvec<-matwhich[modin,]
  #number of parameters
  npar<-sum(parvec)
  #calculates the mean predicted square error
  print(paste("MPSE=",ls.strB$press/(length(yvec)-(npar+1))))
  print(paste("Cp=",z2[modin]))
  print(" Model Statistics ")
  ls.print(ls.strB)
  #predictive analysis 
  print("coeff without int")
  print(ls.strB$beta[-1])
  #plots the predicted y hat against y vec
  plot(xmat[,matwhich[modin,]]%*%ls.strB$beta[-1],yvec,main="Predicted vs Actual",xlab="Predicted",ylab="Actual")
  #calculates the correlation between predicted and actual
  cor<-cor(xmat[,matwhich[modin,]]%*%ls.strB$beta[-1],yvec)
  print(paste("correlation",cor))
  outlist1<-list(correlation=cor)
  outlist1
}

lars.select<-function(xmat,y,int=F,ncheck=10)
{
  #calculates the lars and stores 
  lasso.str<-lars(xmat,y,intercept=int)
  #plots the lasso
  plot(lasso.str)
  #plots the lasso df and lasso cp to see which model has cp about equal to df indicaticing a good model  
  plot(lasso.str$df,lasso.str$Cp,log="y")
  lines(lasso.str$df,lasso.str$df) 
  #plots the cross validated MSE with the fraction of final L1 norm 
  cv.str<-cv.lars(xmat,y,plot.it=T,intercept=int)
  #orders cross validated mse stats and stores the order
  o1<-order(cv.str$cv)
  #stores the index of the smallest cv value
  mindex<-cv.str$index[o1][1]
  #extracts coeffecients from lars 
  beta<-coef(lasso.str)
  #finds the sum of the absolute value of every row in the beta matrix
  index0<-apply(beta,1,sumabs)
  #divide each number by the max index and orders them 
  index0<-index0/max(index0)
  o1<-order(abs(index0-mindex))
  #find minimum index of I1
  I1<-(abs(index0-mindex)==min(abs(index0-mindex)))
  n1<-length(beta[,1])
  beta.out<-beta[I1,]
  #orders using Cp if sumabs is zero 
  if(sum(abs(beta.out))==0){
    v1<-lasso.str$Cp
    o2<-order(v1)
    beta.out<-beta[o1[1:ncheck],]
  }
  Ind.out<-beta.out!=0
  #outputs
  print(beta.out)
  print(Ind.out)
  print(paste("Cp=",lasso.str$Cp[I1]))
  #outlist<-list(beta.out=beta.out,ind.out=Ind.out,Cp=lasso.str$Cp[I1])
  if(int){
    Int.out1<-mean(y)-mean(xmat%*%beta.out[i])
    print(beta.out)
    print(Ind.out)
    print(Ind.out1)
    print(paste("Cp=",lasso.str$Cp[I1]))
  }       
  #plots predicted vs actual 
  plot(xmat%*%beta.out,y,main="Predicted vs Actual",xlab="Predicted",ylab="Actual")
  #calculates correlation
  cor<-cor(xmat%*%beta.out,y)
  print(paste("correlation",cor))
  outlist2<-list(correlation=cor)
  outlist2
}

bestmodelcomp<-function(xmat,y,run.leaps = T,run.lars=T)
{
  if(run.leaps)
  {
    print("running leaps")
    hi<-leaps.then.press(xmat,y)
    corleaps<-hi$correlation
  }
  if(run.lars)
  {
    print("running lars")
    lol<-lars.select(xmat,y)
    corlars<-lol$correlation
  }
  #compares correlations of both models if both are True 
  if(run.leaps&&run.lars)
  {
    if (corleaps<corlars)
    {
      print("Lars method produces the best model")
    }
    else
    {
      print("Leaps method produces the best model")
    }
  }
}

library(ISLR)
library(lars)
library(leaps)

Auto$mpg <- as.numeric(Auto$mpg)

auto.mat<-as.matrix(Auto[,-9])
y.auto<-data.matrix(as.numeric(auto.mat[,1]))
x.auto<-data.matrix(auto.mat[,2:7])
x.auto2<-data.matrix(matrix.2ndorder.make(x.auto))

auto.mat.japan<-data.matrix(auto.mat[auto.mat[,8]==3,])
auto.mat.germany<-data.matrix(auto.mat[auto.mat[,8]==2,])
auto.mat.usa<-data.matrix(auto.mat[auto.mat[,8]==1,])

y.auto.japan<-data.matrix(auto.mat.japan[,1])
x.auto.japan<-data.matrix(auto.mat.japan[,2:7])
x.auto2.japan<-data.matrix(matrix.2ndorder.make(x.auto.japan))

y.auto.germany<-data.matrix(auto.mat.germany[,1])
x.auto.germany<-data.matrix(auto.mat.germany[,2:7])
x.auto2.germany<-data.matrix(matrix.2ndorder.make(x.auto.germany))

y.auto.usa<-data.matrix(auto.mat.usa[,1])
x.auto.usa<-data.matrix(auto.mat.usa[,2:7])
x.auto2.usa<-data.matrix(matrix.2ndorder.make(x.auto.usa))

print("==========================")
print("AUTO (ALL COUNTRIES) MODEL COMPARISON")
print("==========================")
bestmodelcomp(x.auto2,y.auto)

print("==========================")
print("USA MODEL COMPARISON")
print("==========================")
bestmodelcomp(x.auto2.usa,y.auto.usa)

print("==========================")
print("Japan MODEL COMPARISON")
print("==========================")
bestmodelcomp(x.auto2.japan,y.auto.japan)

print("==========================")
print("Germany MODEL COMPARISON")
print("==========================")
bestmodelcomp(x.auto2.germany,y.auto.germany)

#source("/Users/prishabhamre/CS PROJECTS/R/StatLearn/Experiment for A4.R")