%% Demo_RTRL2E
clc,clear

X = double(imread('1.jpg'))./255;
O = X;maxP = max(abs(O(:)));
[n1,n2,n3] = size(X);

CR = 0.5;SR = 1 - CR; 
totalNum = n1*n2*n3;
Omegac =  randperm(totalNum,floor(totalNum*CR)); 
X(Omegac(:)) = 0;
Allindex = 1:totalNum;Omega = setdiff(Allindex,Omegac);

NR = 1;meanmu=0.2;sigma2 = 0.1; 
NRindex = randperm(length(Omega),floor(NR*length(Omega)));
NROmega = Omega(NRindex);
X(NROmega) = X(NROmega) + normrnd(meanmu,sigma2,1,length(NROmega));

K = 500;
tol = 1e-5;
XO = X(Omega);
lambda = 10*max(XO)/sqrt(mean(size(X))*(1-SR));
R = round(0.1*min(n1,n2));
[Xin,Uin,Vin] = RE(X,Omega,Omegac,lambda,R,K,tol);


itmax = 500;
W = ones(n1,n2,n3);
W(Omegac(:)) = 0;
eps = 1e-5;
qta = 0.094;
beta = qta*(1-SR)/(1+SR) *min(n1,n2) * exp(min(n1,n2)/max(n1,n2));
S = R;
etamax = log(50);
lambda = 1200;

mu = 4.5;rho = 0.4;
tic
L = demo_RTRL2E(X,Xin,Uin,Vin,W,etamax,beta,mu,rho,lambda,Omega,S,itmax,eps);
time_RTRL2E = toc;
X_hat_RTRL2E = max(L,0);
X_hat_RTRL2E = min(X_hat_RTRL2E,maxP);
psnr_RTRL2E = PSNR(O,X_hat_RTRL2E,maxP) 

