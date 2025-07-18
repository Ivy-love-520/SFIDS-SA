function [gBestScore,gBest,cg_curve,AK]=PSO(N,Max_iteration,lb,ub,dim,fobj)
%% ��������
% ����λ�ú��ٶȸ��²���
c1=1.5;
c2=1.5; 
maxgen=Max_iteration; % ��������
sizepop=N; % ��Ⱥ��ģ����������
% Ŀ�꺯����Ϣ
fun=fobj;
xub=ub; % ����ȡֵ�Ͻ�
xlb=lb; % ����ȡֵ�½�
vub=[1,1]; % �ٶ��Ͻ�
vlb=[-1,-1]; % �ٶ��½�
%% ��ʼ��
% ����λ�ó�ʼ��
xRange=repmat((xub-xlb),[sizepop,1]);
xLower=repmat(xlb,[sizepop,1]);
pop=rand(sizepop,dim).*xRange+xLower;
% �����ٶȳ�ʼ��
vRange=repmat((vub-vlb),[sizepop,1]);
vLower=repmat(vlb,[sizepop,1]);
V=rand(sizepop,dim).*vRange+vLower;
% �����ʼ����Ŀ�꺯��ֵ����Ӧ�Ⱥ���ֵ��
fitness=ones(sizepop,1);
for k=1:sizepop
    fitness(k)=fun(pop(k,:));
end
%% Ѱ�ҳ�ʼ��ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:); % ȫ�����Ž��λ��
gbest=pop; % ����λ��
fitnessgbest=fitness; % ������Ӧ��ֵ
fitnesszbest=bestfitness; % ȫ��������Ӧ��ֵ
%% ����Ѱ��
yy=ones(sizepop,1); % ���ڱ���ÿ�ε���������Ŀ�꺯��ֵ
for k=1:maxgen
    % ����λ�ú��ٶȸ���
    for m=1:sizepop
        % �����ٶȸ���
        sol=V(m,:)+c1*rand*(gbest(m,:)-pop(m,:))+c2*rand*(zbest-pop(m,:));
        % ȷ�������ٶ�ȡֵ��Χ��Խ��
        ind=find(sol<vlb);
        sol(ind)=vlb(ind);
        ind=find(sol>vub);
        sol(ind)=vub(ind);
        V(m,:)=sol;
        % ����λ�ø���
        sol=pop(m,:)+0.5*V(m,:);
        % ȷ������λ��ȡֵ��Χ��Խ��
        ind=find(sol<xlb);
        sol(ind)=xlb(ind);
        ind=find(sol>xub);
        sol(ind)=xub(ind);
        pop(m,:)=sol;
        % ����������Ӧ��ֵ
        fitness(m)=fun(pop(m,:));
    end
    
    % ���弫ֵ��λ�ú�Ⱥ�弫ֵ��λ�ø���
    for m=1:sizepop
        % ���弫ֵ����λ�ø���
        if fitness(m)<fitnessgbest(m)
            gbest(m,:)=pop(m,:);
            fitnessgbest(m)=fitness(m);
        end
        % Ⱥ�弫ֵ����λ�ø���
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