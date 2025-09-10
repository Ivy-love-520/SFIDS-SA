function [gBestScore,gBest,cg_curve,AK]=PSO(N,Max_iteration,lb,ub,dim,fobj)
%% 参数设置
% 粒子位置和速度更新参数
c1=1.5;
c2=1.5; 
maxgen=Max_iteration; % 迭代次数
sizepop=N; % 种群规模（粒子数）
% 目标函数信息
fun=fobj;
xub=ub; % 参数取值上界
xlb=lb; % 参数取值下界
vub=[1,1]; % 速度上界
vlb=[-1,-1]; % 速度下界
%% 初始化
% 粒子位置初始化
xRange=repmat((xub-xlb),[sizepop,1]);
xLower=repmat(xlb,[sizepop,1]);
pop=rand(sizepop,dim).*xRange+xLower;
% 粒子速度初始化
vRange=repmat((vub-vlb),[sizepop,1]);
vLower=repmat(vlb,[sizepop,1]);
V=rand(sizepop,dim).*vRange+vLower;
% 计算初始化的目标函数值（适应度函数值）
fitness=ones(sizepop,1);
for k=1:sizepop
    fitness(k)=fun(pop(k,:));
end
%% 寻找初始极值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:); % 全局最优解的位置
gbest=pop; % 个体位置
fitnessgbest=fitness; % 个体适应度值
fitnesszbest=bestfitness; % 全局最优适应度值
%% 迭代寻优
yy=ones(sizepop,1); % 用于保存每次迭代的最优目标函数值
for k=1:maxgen
    % 粒子位置和速度更新
    for m=1:sizepop
        % 粒子速度更新
        sol=V(m,:)+c1*rand*(gbest(m,:)-pop(m,:))+c2*rand*(zbest-pop(m,:));
        % 确保粒子速度取值范围不越界
        ind=find(sol<vlb);
        sol(ind)=vlb(ind);
        ind=find(sol>vub);
        sol(ind)=vub(ind);
        V(m,:)=sol;
        % 粒子位置更新
        sol=pop(m,:)+0.5*V(m,:);
        % 确保粒子位置取值范围不越界
        ind=find(sol<xlb);
        sol(ind)=xlb(ind);
        ind=find(sol>xub);
        sol(ind)=xub(ind);
        pop(m,:)=sol;
        % 更新粒子适应度值
        fitness(m)=fun(pop(m,:));
    end
    
    % 个体极值及位置和群体极值及位置更新
    for m=1:sizepop
        % 个体极值及其位置更新
        if fitness(m)<fitnessgbest(m)
            gbest(m,:)=pop(m,:);
            fitnessgbest(m)=fitness(m);
        end
        % 群体极值及其位置更新
        if fitness(m)<fitnesszbest
            zbest=pop(m,:);
            fitnesszbest=fitness(m);
           
         
           
        end
    end
    AK(k,:)=zbest;
  cg_curve(k)= fitnesszbest;
end
gBestScore=fitnesszbest;
gBest=zbest;
end