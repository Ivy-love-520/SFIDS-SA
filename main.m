clear all ;
clc;
close all;
%% 导入数据
fs=5120;%%采样频率，用于频谱图
% 读取纯数值数据
data = load('E01_240731_104328334.txt');

% 转置数据
transposedData = data';

% 保存转置后的数据到新文件
dlmwrite('data.txt', transposedData);
signal=load('data.txt');
% 
% % ssignal1=signal(:,2);
s2=signal(3,:);
% % s3=signal(:,4);
% % s4=signal(:,5);
% % AM1=abs(fft(signal(:,2)));
% % AM2=abs(fft(signal(:,3)));
% % AM3=abs(fft(signal(:,4)));
% 
DATA.data=s2;

% DATA.data=xlsread('data.xlsx');%%读取数据
DATA.tau=0;                    
DATA.dc=0;             
DATA.init=1;           
DATA.tol=1e-9;  
%% 更换适应度函数
fitness = @(x)fitness1(x,DATA);%%综合评价指标1
% fitness = @(x)fitness2(x,DATA);%%包络熵
% fitness = @(x)fitness3(x,DATA);%%包络谱峭度
% fitness = @(x)fitness4(x,DATA);%%幅值谱熵
% fitness = @(x)fitness5(x,DATA);%%模糊熵
% fitness = @(x)fitness6(x,DATA);%%皮尔逊相关系数
% fitness = @(x)fitness7(x,DATA);%% 峭度值
% fitness = @(x)fitness8(x,DATA);%% 样本熵
% fitness = @(x)fitness9(x,DATA);%% 排列熵
% fitness = @(x)fitness10(x,DATA);%% 信息熵

% fitness = @(x)fitness11(x,DATA);%%综合评价指标2
% fitness = @(x)fitness12(x,DATA);%% 多尺度排列熵
% fitness = @(x)fitness13(x,DATA);%% 多尺度样本熵
% fitness = @(x)fitness14(x,DATA);%% 多尺度模糊熵
%% 优化寻参
SearchAgents_no =40;                   % 数量
Max_iteration = 10;                    % 最大迭代次数
dim = 2;                               % 优化参数个数
lb = [1000, 8];                 % 参数取值下界(a,k)
ub = [4000, 15];                 % 参数取值上界
%% 参数寻优
[Best_score,Best_pos,Convergence_curve,ak]=PSO(SearchAgents_no,Max_iteration,lb ,ub,dim,fitness);
%% 显示最优参数并代入VMD
disp(['最优K值为：',num2str(round(Best_pos(1,2)))])
disp(['最优alpha值为：',num2str(round(Best_pos(1,1)))])
disp(['最优指标为：',num2str(Best_score)])
[u, u_hat, omega] = VMD(DATA.data, round(Best_pos(1,1)), DATA.tau, round(Best_pos(1,2)), DATA.dc, DATA.init, DATA.tol);
%% 可视化
%% 适应度曲线
figure
plot(Convergence_curve,'linewidth',1.5);
title('收敛曲线')
xlabel('迭代次数')
ylabel('适应度值')
xlim([1,Max_iteration]);
grid on
saveas(gcf,['结果图/1',],'fig');
[m,n]=size(u);
%% 分解结果可视化
figure
subplot(m+1,1,1);
plot(DATA.data,'k');grid on;
title('原始图像');
for i = 1:m
    subplot(m+1,1,i+1);
    plot(u(i,:),'k');
    title(['PSO-VMD分解信号',num2str(i),':']);
end
saveas(gcf,['结果图/2',],'fig');
set(gcf,'position',[50 50 650 650])%%调整画布大小
%% 频谱图可视化
figure
for i=1:m
subplot(m,1,i)
%% FFT 变换
[cc,y_f]=hua_fft(u(i,:),fs,1);
plot(y_f,cc,'b','LineWIdth',1.5);
% hua_fft(u(i,:),fs,1)
ylabel(['FFT of IMF',num2str(i)]);
axis tight
end
saveas(gcf,['结果图/3',],'fig');
%% 最优参数变化曲线
figure
plot(round(ak(:,1)),'k-s ','linewidth',1);
title('惩罚因子优化过程曲线')
xlabel('迭代次数')
ylabel('惩罚因子')
xlim([1,Max_iteration]);
ylim([lb(:,1),ub(:,1)]);
saveas(gcf,['结果图/4',],'fig');
figure
plot(round(ak(:,2)),'k-d ','linewidth',1);
title('分解模态数优化过程曲线')
xlabel('迭代次数')
ylabel('分解模态数')
xlim([1,Max_iteration]);
ylim([lb(:,2),ub(:,2)]);
saveas(gcf,['结果图/5',],'fig');
%% 分解结果3D视图
u1=u(1:round(Best_pos(1,2)),:);
plot3imf(u1);%%为简便以函数调用，函数在plot3imf文件
title('分解结果3D视图')
saveas(gcf,['结果图/6',],'fig');

save 结果图/res
