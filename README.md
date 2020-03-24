# Logistic regression and Poisson regression
```
YY<-c(0 ,0,1,1,1,2,0,0,1,2,
      0 ,1,1,0,2,0,0,3,0,2,
      1 ,0,0,0,0,0,4,0,0,1,
      12,1,0,2,2,4,0,1,0,1)

XX<-cbind(1,rep(c(0,0.1,0.2,0.3),each=10))

nn<-c(15, 3, 9,12,13,13,16,11,11, 8,
      6,14,12,10,14,12,14,14,10,12,
      12,12,11,13,12,14,15,14,12, 6,
      12,12,13, 8,12,13,13,13,12, 9)

beta<-rep(0,2)
for(iter in 1:100){
  mu<-exp(XX%*%beta+log(nn))
  vv<-mu
  Delta<-1/mu
  WW<-1/vv/(Delta)^2
  UU<-t(XX)%*%diag(c(WW*Delta))%*%(YY-mu)
  II<-t(XX)%*%diag(c(WW))%*%XX
  beta<-beta+solve(II,UU)
}
beta

beta<-rep(0,2)

for(iter in 1:100){
  mu<-exp(XX%*%beta+log(nn))
  vv<-mu
  Delta<-1/mu
  WW<-1/vv/(Delta)^2
  zz<-XX%*%beta+c(Delta)*c(YY-mu)
  beta<-solve(t(XX)%*%diag(c(WW))%*%XX,t(XX)%*%diag(c(WW))%*%zz)
}
beta

glm(YY~-1+XX+offset(log(nn)),family=poisson())

##########

pp<-YY/nn
temp<-1/nn

glm(cbind(YY,nn-YY)~-1+XX,family=binomial())
glm(pp~-1+XX,weights=nn,family=binomial())

beta<-rep(0,2)
for(iter in 1:100){
  mu<-exp(XX%*%beta)/(1+exp(XX%*%beta))
  vv<-temp*mu*(1-mu)
  Delta<-1/(mu*(1-mu))
  WW<-1/vv/(Delta)^2
  UU<-t(XX)%*%diag(c(WW*Delta))%*%(pp-mu)
  II<-t(XX)%*%diag(c(WW))%*%XX
  beta<-beta+solve(II,UU)
}
beta

beta<-rep(0,2)
for(iter in 1:100){
  mu<-exp(XX%*%beta)/(1+exp(XX%*%beta))
  vv<-temp*mu*(1-mu)
  Delta<-1/(mu*(1-mu))
  WW<-1/vv/(Delta)^2
  zz<-XX%*%beta+c(Delta)*c(pp-mu)
  beta<-solve(t(XX)%*%diag(c(WW))%*%XX,t(XX)%*%diag(c(WW))%*%zz)
}
beta
```
