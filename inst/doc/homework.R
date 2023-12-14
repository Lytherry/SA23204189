## ----results='hide'-----------------------------------------------------------
library(alr4)
library(ggplot2)

## -----------------------------------------------------------------------------
str(brains)#查看数据集brains的数据类型、变量名称和观测值

## -----------------------------------------------------------------------------
summary(brains)#显示脑重、体重的均值、分位数

## -----------------------------------------------------------------------------
model<-lm(log10(BrainWt)~log10(BodyWt),data=brains)
summary(model)

## -----------------------------------------------------------------------------
ggplot(data=brains,aes(y=BrainWt,x=BodyWt))+
  scale_x_log10()+scale_y_log10()+                                   #将坐标尺度修改为每格10倍
  geom_point(color="red")+                                           #画点
  stat_smooth(method=lm,formula=y~x,color="blue")+                   #画回归直线
  ylab("脑重")+xlab("体重")+                                         #修改坐标名称
  ggtitle("哺乳动物脑重和体重的拟合图像")+                           #加标题
  theme(plot.title = element_text(hjust=0.5))+
  annotate("text",label=expression(R^2==0.9208),x=1000,y=50)+        #加图例
  annotate("text",label=expression(log[10](y)==0.92713+0.75169%.%log[10](x)),x=500,y=10)

## -----------------------------------------------------------------------------
BodyWt_min<-min(brains$BodyWt)
BodyWt_max<-max(brains$BodyWt)
predicted<-data.frame(BodyWt=seq(BodyWt_min,BodyWt_max,length.out=100)) #创建需要预测的自变量数值
predicted$BrainWt <- 10^(predict(model,predicted)) #预测数值
head(predicted)

## -----------------------------------------------------------------------------

plot(model)

## -----------------------------------------------------------------------------
u<-runif(1000)      #生成1000个标准柯西变量
x<-tan(pi*u-pi/2)
head(x)             #显示头几个柯西变量的值

## -----------------------------------------------------------------------------
hist(x[-10<x&x<10],prob=T,xlim=c(-10,10),ylim=c(0,0.35),breaks=30,main="样本直方图与密度函数对比图",col="yellow")
curve(1/pi*1/(1+x^2),col="red",add=T,lwd=2)

## -----------------------------------------------------------------------------
my.sample<-function(x,size,prob=rep(1/length(x),length(x))){  
  #x代表样本空间，size代表样本大小，prob代表概率(缺省值情况下为等概率)
  
  stopifnot(length(x)==length(prob))
  #x的长度应当和prob的长度相等
  
  cp<-cumsum(prob)
  U<-runif(size)              #生成均匀分布样本
  R<-x[findInterval(U,cp)+1]  #生成所求样本
  return(R)
}

## -----------------------------------------------------------------------------
set.seed(1)
x<-c(1,2,"a",TRUE)
sample(x,30,replace=TRUE,prob=c(0.1,0.2,0.3,0.4))

## -----------------------------------------------------------------------------
set.seed(2)
my.sample(x,30,prob=c(0.1,0.2,0.3,0.4))

## -----------------------------------------------------------------------------
pi_est<-function(rho,n,K){   #rho=l/d,n为针数，K为重复模拟的次数
  l<-rho
  d<-1
  pi_hat<-numeric(K)
  
  for(i in 1:K){
  X<-runif(n,0,d/2)  #针心到平行线的最近距离
  Y<-runif(n,0,pi/2) #针与平行线的夹角
  pi_hat[i]<-2*l/d/mean(0.5*l*sin(Y)>X)
  }
  
  pi_hat
}

## -----------------------------------------------------------------------------
n<-10^6
K<-100

set.seed(1)
pi_hat1<-pi_est(0.5,n,K)
pi_hat2<-pi_est(0.8,n,K)
pi_hat3<-pi_est(1,n,K)
print(paste("rho=0.5时的pi估计值为",mean(pi_hat1),"方差为",var(pi_hat1)),sep="")
print(paste("rho=0.8时的pi估计值为",mean(pi_hat2),"方差为",var(pi_hat2)),sep="")
print(paste("rho=1时的pi估计值为",mean(pi_hat3),"方差为",var(pi_hat3)),sep="")

## -----------------------------------------------------------------------------
n<-10000  #每次取10000个数估计积分值
M<-100    #将积分值估计100次
Y1<-numeric(M)
Y2<-numeric(M)

set.seed(8)
for(i in 1:M){
X<-runif(n)

#简单蒙特卡罗方法
y1<-exp(X)
Y1[i]<-mean(y1)

#对偶变量法
y2<-c(exp(X[1:n/2]),exp(1-X[1:n/2]))
Y2[i]<-mean(y2)
}

## -----------------------------------------------------------------------------
mean(Y1) #简单蒙特卡罗方法估计量均值
var(Y1)  #简单蒙特卡罗方法估计量方差
mean(Y2) #对偶变量法估计量均值
var(Y2)  #对偶变量法估计量方差

print(paste("该实验中，对偶变量法将方差缩减了",round(100*(1-var(Y2)/var(Y1)),2),"%",sep=""))

## -----------------------------------------------------------------------------
m <- 1e6
est <- sd <- numeric(2)
g <- function(x) {
  x^2 * exp(-x^2/2) / sqrt(2*pi) * (x > 1) 
}

set.seed(1)

#构造f1
f1<-function(x){
  exp(-x^2/2) / (1-pnorm(1)) /sqrt(2*pi)
}

#生成m个服从f1的随机变量
x<-numeric(m)
item<-1
while(item<=m){
  x[item]<-rnorm(1)
  if(x[item]>1)
    item<-item+1
}

fg<-g(x)/f1(x)
est[1]<-mean(fg)
sd[1]<-sd(fg)

#构造f2
f2<-function(x){
  1/(x^2)
}

#生成m个服从f2的随机变量(这里运用了逆函数法)
u<-runif(m)
x<-1/(1-u)

fg<-g(x)/f2(x)
est[2]<-mean(fg)
sd[2]<-sd(fg)

#估计值的均值
est
#估计值的标准差
sd

## ----fig.width=10-------------------------------------------------------------
    x <- seq(1.01, 5, .01)
    w<-2
    gs<-c("g","f1","f2")

    #for color change lty to col
    #figure (a)
    plot(x, g(x), type = "l", ylab = "",
         ylim = c(0,2), lwd = w, col=1,main='(A)')
    lines(x, f1(x), lty = 2, lwd = w,col=2)
    lines(x, f2(x), lty = 3, lwd = w,col=3)
    legend("topright", legend = gs,
           lty = 1:3, lwd = w, inset = 0.02,col=1:3)

    #figure (b)
    plot(x, rep(1,length(x)), type = "l", ylab = "",
        ylim = c(0,3.2), lty = 1, lwd = w,col=1,main='(B)')
    lines(x, g(x)/f1(x), lty = 2, lwd = w,col=2)
    lines(x, g(x)/f2(x), lty = 3, lwd = w,col=3)
    legend("topright", legend = gs,
           lty = 1:3, lwd = w, inset = 0.02,col=1:3)

## -----------------------------------------------------------------------------
M <- 10000; k <- 5 # 分成5个子区间
r <- M/k #每个区间取的随机数个数
N <- 50 #实验次数
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)

set.seed(2)

for (i in 1:N) {
  #使用例5.10中的方法
  u<-runif(M)
  x<- -log(1-u*(1-exp(-1)))
  fg<-g(x)/exp(-x)*(1-exp(-1))
  est[i,1]<-mean(fg)
  
  #使用分层抽样法
  for(j in 1:k){
    u<-runif(r,(j-1)/k,j/k)
    x<- -log(1-u*(1-exp(-1)))
    fg[((j-1)*r+1):(j*r)]<-g(x)/exp(-x)*(1-exp(-1))
  }
  est[i,2]<-mean(fg)
}

#两种方法的均值对比
round(apply(est,2,mean),5)
#两种方法的标准差对比
round(apply(est,2,sd),5)

## -----------------------------------------------------------------------------
m<-100000 #重复实验次数
n<-20   #每次实验取20个样本
df<-2   #自由度为2的卡方分布
alpha<-0.05
d1<-numeric(m) #置信下限
d2<-numeric(m) #置信上限

set.seed(3)

for(i in 1:m){
  x<-rchisq(n,df)
  d1[i]<-mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  d2[i]<-mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
}

mu0<-2
mean(d1<mu0&d2>mu0) #mu0落在置信区间的概率

## -----------------------------------------------------------------------------
m<-100000 #重复实验次数
n<-20   #每次实验取20个样本
df<-1   #自由度为1的卡方分布
alpha<-0.05
d1<-numeric(m) #置信下限
d2<-numeric(m) #置信上限

set.seed(4)

for(i in 1:m){
  x<-rchisq(n,df)
  d1[i]<-mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  d2[i]<-mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
}

mu0<-1
1-mean(d1<mu0&d2>mu0) #mu0落在置信区间外的概率

## -----------------------------------------------------------------------------
m<-100000 #重复实验次数
n<-20   #每次实验取20个样本
alpha<-0.05
d1<-numeric(m) #置信下限
d2<-numeric(m) #置信上限

set.seed(5)

for(i in 1:m){
  x<-runif(n,0,2)
  d1[i]<-mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  d2[i]<-mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
}

mu0<-1
1-mean(d1<mu0&d2>mu0) #mu0落在置信区间外的概率

## -----------------------------------------------------------------------------
m<-100000 #重复实验次数
n<-20   #每次实验取20个样本
alpha<-0.05
d1<-numeric(m) #置信下限
d2<-numeric(m) #置信上限

set.seed(6)

for(i in 1:m){
  x<-rexp(n,1)
  d1[i]<-mean(x)-sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
  d2[i]<-mean(x)+sd(x)*qt(1-alpha/2,n-1)/sqrt(n)
}

mu0<-1
1-mean(d1<mu0&d2>mu0) #mu0落在置信区间外的概率

## -----------------------------------------------------------------------------
options(max.print=3000)

m<-1000
m1<-950
m2<-50

set.seed(1)
p<-numeric(m)                 #生成1000个p值
p[1:m1]<-runif(m1,0,1)        #前950个变量服从(0,1)均匀分布
p[(m1+1):m]<-rbeta(m2,0.1,1)  #后50个变量服从beta(0.1,1)

p.sorted<-sort(p)
p.ordered<-order(p)

length(p)
head(p)
head(p.sorted)
head(p.ordered)

## -----------------------------------------------------------------------------
p.adj1 <- p.sorted/(1:m)*m
p.adj2 <- p.adjust(p.sorted,method='BH')
p.adj3 <- p.adjust(p.sorted,method='bonferroni')
rbind(p.adj1,p.adj2,p.adj3)[,1:50]

## -----------------------------------------------------------------------------
alpha<-0.1
p.refuse2<-(p.adj2>0.1) #TRUE为接受原假设，FALSE为拒绝原假设
p.refuse3<-(p.adj3>0.1)

## -----------------------------------------------------------------------------
A<-matrix(nrow=2,ncol=3)
colnames(A)<-c("FWER","FDR","TPR")
rownames(A)<-c("Bonf","B-H")
A[1]<-1-mean(p.refuse3)
A[2]<-1-mean(p.refuse2)
A[3]<-length(which(p.refuse3==FALSE&p.ordered<=950))/length(which(p.refuse3==FALSE))
A[4]<-length(which(p.refuse2==FALSE&p.ordered<=950))/length(which(p.refuse2==FALSE))
A[5]<-length(which(p.refuse3==FALSE&p.ordered>950))/m2
A[6]<-length(which(p.refuse2==FALSE&p.ordered>950))/m2
A

## -----------------------------------------------------------------------------
bootstrap<-function(lambda,n,B,m){ #B为每次模拟生成bootstrap变量的个数，m为模拟次数
  bias<-numeric(m)    #bootstrap均值偏差
  se.boot<-numeric(m) #bootstrap标准差
  se.samp<-numeric(m) #样本标准差
  
  for(i in 1:m){
    x<-rexp(n,lambda)
    theta<-1/mean(x)
    thetastar<-numeric(B)
    for(b in 1:B){
      xstar<-sample(x,replace=TRUE)
      thetastar[b]<-1/mean(xstar)
    }
    bias[i]<-mean(thetastar)-theta
    se.boot[i]<-sd(thetastar)
    se.samp[i]<-sd(x)/sqrt(length(x))
  }
  
  data.frame(bias,se.boot,se.samp)
}

## -----------------------------------------------------------------------------
lambda<-2
n1<-5
n2<-10
n3<-20
B<-1000
m<-1000

set.seed(2)
bootstrap1<-bootstrap(lambda,n1,B,m)
bootstrap2<-bootstrap(lambda,n2,B,m)
bootstrap3<-bootstrap(lambda,n3,B,m)

## -----------------------------------------------------------------------------
mean(bootstrap1$bias) #理论值为0.5
mean(bootstrap2$bias) #理论值为0.222
mean(bootstrap3$bias) #理论值为0.105

## -----------------------------------------------------------------------------
mean(bootstrap1$se.boot) #理论值为1.443
mean(bootstrap2$se.boot) #理论值为0.786
mean(bootstrap3$se.boot) #理论值为0.496

## -----------------------------------------------------------------------------
library("bootstrap")

print(cor(law$LSAT,law$GPA))
print(cor(law82$LSAT,law82$GPA))

## -----------------------------------------------------------------------------
B<-200
alpha<-0.05
n<-nrow(law)
thetastar<-numeric(B)     #theta的bootstrap估计
se.thetastar<-numeric(B)  #se.theta的bootstrap估计

theta<-cor(law$LSAT,law$GPA)
se.theta<-0.115

set.seed(3)

for(b in 1:B){
  #theta的bootstrap估计
  i<-sample(1:n,size=n,replace=TRUE)
  LSAT<-law$LSAT[i]
  GPA<-law$GPA[i]
  thetastar[b]<-cor(LSAT,GPA)
  
  #se.theta的bootstrap估计
  tempstar<-numeric(B)
  for(c in 1:B){
    j<-sample(1:n,size=n,replace=TRUE)
    tempstar[c]<-cor(law$LSAT[j],law$GPA[j])
  }
  se.thetastar[b]<-sd(tempstar)
}

tstar<-unname(quantile((thetastar-theta)/se.thetastar,1-alpha/2))  #t的1-alpha/2分位数,alpha=0.05
d1<-theta-tstar*se.theta  #bootstrap置信下限估计
d2<-theta+tstar*se.theta  #bootstrap置信上限估计

tstar
d1
d2

## -----------------------------------------------------------------------------
library(boot)
m<-1e3
boot.mean<-function(x,i) mean(x[i])
sample<-c(3,5,7,18,43,85,91,98,100,130,230,487)
ci.norm<-matrix(NA,m,2)   #标准正态法
ci.basic<-matrix(NA,m,2)  #基本法
ci.perc<-matrix(NA,m,2)   #百分位数法
ci.bca<-matrix(NA,m,2)    #BCa方法
set.seed(1)

for(i in 1:m){
  R<-sample(sample,replace=TRUE)
  de<-boot(data=R,statistic=boot.mean,R=999)
  ci<-boot.ci(de,type=c("norm","basic","perc","bca"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
  ci.bca[i,]<-ci$bca[4:5]
}

cat('norm: [',mean(ci.norm[,1]),',',mean(ci.norm[,2]),']','\n',
'basic: [',mean(ci.basic[,1]),',',mean(ci.basic[,2]),']','\n',
'perc: [',mean(ci.perc[,1]),',',mean(ci.perc[,2]),']','\n',
'BCa: [',mean(ci.bca[,1]),',',mean(ci.bca[,2]),']',sep="")

## -----------------------------------------------------------------------------
library("bootstrap")
n<-nrow(scor)

cov.hat<-cov(scor)
eigen.hat<-sort(eigen(cov.hat)$values,decreasing = TRUE)
theta.hat<-eigen.hat[1]/sum(eigen.hat)

theta.jack<-numeric(n)
for(i in 1:n){
  cov.jack<-cov(scor[-i,])
  eigen.jack<-sort(eigen(cov.jack)$values,decreasing = TRUE)
  theta.jack[i]<-eigen.jack[1]/sum(eigen.jack)
}

bias.jack<-(n-1)*(mean(theta.jack)-theta.hat)
se.jack<-sqrt((n-1)/n*sum((theta.jack-theta.hat)^2))
c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)

n<-length(magnetic)
e1<-e2<-e3<-e4<-matrix(data=0,nrow=2,ncol=n*(n-1)/2)   #每种方法的偏差，用四个矩阵存储
yhat1<-yhat2<-logyhat3<-logyhat4<-numeric(2)
col<-0                                          #将数据输入矩阵的第col列

for(i in 1:(n-1)){
  for(j in (i+1):n){
    col<-col+1
    y<-magnetic[c(-i,-j)]
    x<-chemical[c(-i,-j)]
    
    #线性
    J1<-lm(y~x)
    yhat1[1]<-J1$coef[1]+J1$coef[2]*chemical[i]
    yhat1[2]<-J1$coef[1]+J1$coef[2]*chemical[j]
    e1[1,col]<-magnetic[i]-yhat1[1]
    e1[2,col]<-magnetic[j]-yhat1[2]
    
    #二次
    J2<-lm(y~x+I(x^2))
    yhat2[1]<-J2$coef[1]+J2$coef[2]*chemical[i]+J2$coef[3]*chemical[i]^2
    yhat2[2]<-J2$coef[1]+J2$coef[2]*chemical[j]+J2$coef[3]*chemical[j]^2
    e2[1,col]<-magnetic[i]-yhat2[1]
    e2[2,col]<-magnetic[j]-yhat2[2]
    
    #指数
    J3<-lm(log(y)~x)
    logyhat3[1]<-J3$coef[1]+J3$coef[2]*chemical[i]
    logyhat3[2]<-J3$coef[1]+J3$coef[2]*chemical[j]
    yhat3<-exp(logyhat3)
    e3[1,col]<-magnetic[i]-yhat3[1]
    e3[2,col]<-magnetic[j]-yhat3[2]
    
    #双对数
    J4<-lm(log(y)~log(x))
    logyhat4[1]<-J4$coef[1]+J4$coef[2]*log(chemical[i])
    logyhat4[2]<-J4$coef[1]+J4$coef[2]*log(chemical[j])
    yhat4<-exp(logyhat4)
    e4[1,col]<-magnetic[i]-yhat4[1]
    e4[2,col]<-magnetic[j]-yhat4[2]
  }
}

c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))

## -----------------------------------------------------------------------------
y<-magnetic
x<-chemical
lm(y~x+I(x^2))

## -----------------------------------------------------------------------------
#计算积分二次方距离
cvm<-function(X,Y){
  n<-length(X)
  m<-length(Y)
  
  X<-sort(X)
  Y<-sort(Y)
  Fn<-function(x) sum(X<=x)/n
  Gm<-function(x) sum(Y<=x)/m
  
  C<-m*n/(m+n)^2
  sigma1<-sum((Fn(X)-Gm(X))^2)
  sigma2<-sum((Fn(Y)-Gm(Y))^2)
  
  return(C*(sigma1+sigma2))
}

#读入数据
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

#进行置换检验
set.seed(1)
options( warn = - 1)
R <- 999           #重复次数
z <- c(x, y)       #原始数据
K <- 1:26
reps <- numeric(R) #存储重复变量
t0 <- cvm(x,y)
for (i in 1:R) {
k <- sample(K, size = 14, replace = FALSE)
x1 <- z[k]
y1 <- z[-k]
reps[i] <- cvm(x1,y1)
}

p <- mean(c(t0, reps) >= t0)
p

## -----------------------------------------------------------------------------
counttest <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy))) #返回端点个数最大值
}

## -----------------------------------------------------------------------------
permutation.counttest<-function(x,y){
R <- 999           #重复次数
z <- c(x, y)       #原始数据
n<-length(x)
m<-length(y)
K<-1:(m+n)

reps <- numeric(R) #存储重复变量
t0 <- counttest(x,y)
for (i in 1:R) {
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]
y1 <- z[-k]
reps[i] <- counttest(x1,y1)
}

p <- mean(c(t0, reps) >= t0)
p
}

## -----------------------------------------------------------------------------
set.seed(3)
x1<-rnorm(50,0,1)
y1<-rnorm(100,0,1)
permutation.counttest(x1,y1)

## -----------------------------------------------------------------------------
set.seed(4)
x2<-rnorm(50,0,1)
y2<-rnorm(100,0,1.5)
permutation.counttest(x2,y2)

## -----------------------------------------------------------------------------
alpha.determination<-function(N,b1,b2,b3,f0){
  X1<-rpois(N,lambda=1)
  X2<-rexp(N,rate=1)
  X3<-rbinom(N,1,prob=0.5)
  
  g<-function(alpha){
    temp<-exp(-alpha-b1*X1-b2*X2-b3*X3)
    q<-1/(1+temp)
    mean(q)-f0
  }
  
  solution<-uniroot(g,c(-100,0))
  unlist(solution)[1]
}

## -----------------------------------------------------------------------------
N<-10^6
b1<-0
b2<-1
b3<--1
f0<-c(0.1,0.01,0.001,0.0001)
n<-length(f0)
a.root<-numeric(n)

set.seed(1)
for(i in 1:n){
  a.root[i]<-alpha.determination(N,b1,b2,b3,f0[i])
}

a.root

## -----------------------------------------------------------------------------
plot(-log(f0),a.root)

## -----------------------------------------------------------------------------
f<-function(x){
  1/2*exp(-abs(x))
}

## -----------------------------------------------------------------------------
mh<-function(sd){
m<-10000
x<-numeric(m)

x[1]<-rnorm(1)
k<-0
u<-runif(m)

for(i in 2:m){
  xt<-x[i-1]
  y<-rnorm(1,mean=xt,sd=sd)
  num<-f(y)*dnorm(xt,mean=y,sd=sd)
  den<-f(xt)*dnorm(y,mean=xt,sd=sd)
  if(u[i]<=num/den){
    x[i]<-y
  }
  else{
    x[i]<-xt
    k<-k+1    #被拒绝的y的个数
  }
}

out<-list(rate=1-k/m,chain=x)
return(out)
}

## -----------------------------------------------------------------------------
set.seed(2)
mh1<-mh(0.5)
mh2<-mh(1)
mh3<-mh(2)

x1<-mh1$chain
x2<-mh2$chain
x3<-mh3$chain

index<-6000:6500

y1<-x1[index]
y2<-x2[index]
y3<-x3[index]

plot(index,y1,type="l",ylab="x")
plot(index,y2,type="l",ylab="x")
plot(index,y3,type="l",ylab="x")

c(mh1$rate,mh2$rate,mh3$rate)

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 5000                #链条长度
burn <- 1000             #链条要去掉前1000个数
X <- matrix(0, N, 2)     #链条（二维变量）

rho <- 0.9 #相关系数
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

#生成链条
set.seed(3)
X[1, ] <- c(mu1, mu2) #初始化
for (i in 2:N) {
x2 <- X[i-1, 2]
m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
X[i, 1] <- rnorm(1, m1, s1)
x1 <- X[i, 1]
m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
X[i, 2] <- rnorm(1, m2, s2)
}

b <- burn + 1
x <- X[b:N, ]
cor(x)

## -----------------------------------------------------------------------------
plot(x, main="", cex=.5, xlab=bquote(X[1]),
ylab=bquote(X[2]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
x_data.frame<-data.frame(x1=x[,1],x2=x[,2])
model<-lm(x2~x1,data=x_data.frame)
summary(model)

## -----------------------------------------------------------------------------
plot(model)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

## -----------------------------------------------------------------------------
f <- function(x, sigma) {
if (any(x < 0)) return (0)
stopifnot(sigma > 0)
return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

chisq.chain<-function(sigma,m){ #sigma=4,m=10000
xt <- x[i-1]
y <- rchisq(1, df = xt)

x <- numeric(m)
x[1] <- rchisq(1, df=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rchisq(1, df = xt)
num <- f(y, sigma) * dchisq(xt, df = y)
den <- f(xt, sigma) * dchisq(y, df = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1 #y is rejected
}
}

return(x)
}

## -----------------------------------------------------------------------------
sigma <- 1 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length

#generate the chains
set.seed(6)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- chisq.chain(sigma, n)

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)

## -----------------------------------------------------------------------------
f<-function(x){
  u<-c(11,8,27,13,16,0,23,10,24,2)
  v<-u+1

  sum((v*exp(-v*x)-u*exp(-u*x))/(exp(-u*x)-exp(-v*x)))
}

## -----------------------------------------------------------------------------
uniroot(f,interval = c(0.01,1),tol=10^-5)

## -----------------------------------------------------------------------------
est<-10^-5                      #误差
n<-10
lambda0<-lambda1<-0.5
temp<-0
iter<-0                         #计数器

while(abs(temp-lambda1)>est){
  lambda1<-n/(n/lambda0-f(lambda0))
  temp<-lambda0
  lambda0<-lambda1
  iter<-iter+1
}

lambda0
iter

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

## -----------------------------------------------------------------------------
#enter the payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
B <- A+2
s <- solve.game(B)

## -----------------------------------------------------------------------------
round(cbind(s$x, s$y), 7)

## -----------------------------------------------------------------------------
x<-c(1,2,3)
y<-matrix(c(1,2,3,4),nrow=2)
dim(x)
dim(y)

## -----------------------------------------------------------------------------
x<-matrix(c(1,2,3,4),nrow=2)
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
x<-data.frame(1:4,letters[1:4])
x<-as.matrix(x)
x
typeof(x)
is.data.frame(x)
is.matrix(x)

## -----------------------------------------------------------------------------
test=data.frame(col1=character(0),col2=numeric(0),col3=logical(0))
str(test)
test

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)  #获取一个长为2的向量，储存x的最小值和最大值
(x - rng[1]) / (rng[2] - rng[1])  #将x的每个元素映射到[0,1]上
}

## -----------------------------------------------------------------------------
data<-data.frame(1:5,c(F,F,F,T,T),2:6)  #构造一个数据框

## -----------------------------------------------------------------------------
#将scale01运用到数据框每一列
lapply(data,scale01)

## -----------------------------------------------------------------------------
#将scale01运用到数据框每一numeric列
lapply(data[lapply(data,is.numeric)==T],scale01)

## -----------------------------------------------------------------------------
#数据框全是numeric型元素
data1<-data.frame(1:5,2:6)  #构造一个纯数值数据框
vapply(data1,sd,numeric(1))

## -----------------------------------------------------------------------------
#数据框是混合型元素
data2<-data.frame(1:5,letters[1:5],2:6)  #构造一个混合型数据框
vapply(data2[vapply(data,is.numeric,logical(1))==T],sd,numeric(1))

## -----------------------------------------------------------------------------
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
a<-2
b<-2
n<-100
x<-1     #x的初始值
y<-0.5   #y的初始值

###### 生成链 #####
set.seed(1)
X[1, ] <- c(x, y) #初始化
for (i in 2:N) {
y <- X[i-1, 2]
X[i, 1] <- rbinom(1,n,y)
x <- X[i, 1]
X[i, 2] <- rbeta(1,x+a,n-x+b)
}

d <- burn + 1
x0 <- X[d:N, ]

plot(x0)

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp)
#  
#  sourceCpp(code='
#  #include <Rcpp.h>
#  #include <iostream>
#  #include <ctime>
#  #include <random>
#  using namespace Rcpp;
#  
#  //[[Rcpp::export]]
#  long cls_random::randomBinomial(long N,double probability){
#      long rnd = 0;
#      for (long i=0;i<N;i++){
#          double pV = (double)rand()/(double)RAND_MAX;
#          if (pV<probability){
#              rnd++;
#          }
#      }
#    return rnd;
#  }
#  
#  double cls_random::randomBeta(double alpha,double beta){
#      double u, v;
#      double x, y;
#      do
#      {
#          u=cls_random::randomUniform();
#          v=cls_random::randomUniform();
#          x=pow(u,1/alpha);
#          y=pow(v,1/beta);
#      } while (x+y>1);
#      return x/(x+y);
#  }
#  
#    NumericVector iters(int a,int b,int n,int x,int y){
#      int N=5000;
#      int X[N][2];
#      int x1;
#      int y1;
#  
#      X[0][0]=x;
#      X[0][1]=y;
#  
#      for(int i=1;i<N;i++) {
#      y1=X[i-1][1];
#      X[i][0]=cls_random::randomBinomial(n,y);
#      x1=X[i][0];
#      X[i][1]=cls_random::randomBeta(x+a,n-x+b);
#      }
#  
#      return X;
#    }
#  ')
#  
#  X<-iters(2,2,100,1,0.5)
#  x0 <- X[1001:10000, ]
#  plot(x0)

