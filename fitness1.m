 %%�ۺ�����ָ��
function fitness1 = fical(Positions,DATA);
X=DATA.data;
alpha = round(Positions(1, 1));
K = round(Positions(1, 2));
tau = DATA.tau;                    
DC = DATA.dc;             
init = DATA.init;           
tol =DATA.tol;     
[u, u_hat, omega] = VMD(X, alpha, tau, K, DC, init, tol);

 %% ��ȡ�ź�����ָ��
    for ii=1:K
         m=2;%ά��һ��ȡ2
    r=0.2*std(u(ii,:));
    tau=1;%%�Ƿ��ز�����1��2�ǣ���ѡ1��
    feature(ii) = sampleEntropy(u(ii,:),m,r,tau);%%������
    end
     Y=sum(u,1);
     pear=corr(X',Y','type','pearson');%%Ƥ��ѷϵ��
    [d,~]=size(omega);%�ۺϴ���
    D=log10(d);%��Ծۺϴ���
    fitness1= min((feature/pear)*D);

    end