 %%综合评价指标
function fitness1 = fical(Positions,DATA);
X=DATA.data;
alpha = round(Positions(1, 1));
K = round(Positions(1, 2));
tau = DATA.tau;                    
DC = DATA.dc;             
init = DATA.init;           
tol =DATA.tol;     
[u, u_hat, omega] = VMD(X, alpha, tau, K, DC, init, tol);

 %% 提取信号评价指标
    for ii=1:K
         m=2;%维数一般取2
    r=0.2*std(u(ii,:));
    tau=1;%%是否重采样，1否2是，常选1；
    feature(ii) = sampleEntropy(u(ii,:),m,r,tau);%%样本熵
    end
     Y=sum(u,1);
     pear=corr(X',Y','type','pearson');%%皮尔逊系数
    [d,~]=size(omega);%聚合代数
    D=log10(d);%相对聚合代数
    fitness1= min((feature/pear)*D);

    end