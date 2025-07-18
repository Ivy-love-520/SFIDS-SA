clear all ;
clc;
close all;
%% ��������
fs=5120;%%����Ƶ�ʣ�����Ƶ��ͼ
% ��ȡ����ֵ����
data = load('E01_240731_104328334.txt');

% ת������
transposedData = data';

% ����ת�ú�����ݵ����ļ�
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

% DATA.data=xlsread('data.xlsx');%%��ȡ����
DATA.tau=0;                    
DATA.dc=0;             
DATA.init=1;           
DATA.tol=1e-9;  
%% ������Ӧ�Ⱥ���
fitness = @(x)fitness1(x,DATA);%%�ۺ�����ָ��1
% fitness = @(x)fitness2(x,DATA);%%������
% fitness = @(x)fitness3(x,DATA);%%�������Ͷ�
% fitness = @(x)fitness4(x,DATA);%%��ֵ����
% fitness = @(x)fitness5(x,DATA);%%ģ����
% fitness = @(x)fitness6(x,DATA);%%Ƥ��ѷ���ϵ��
% fitness = @(x)fitness7(x,DATA);%% �Ͷ�ֵ
% fitness = @(x)fitness8(x,DATA);%% ������
% fitness = @(x)fitness9(x,DATA);%% ������
% fitness = @(x)fitness10(x,DATA);%% ��Ϣ��

% fitness = @(x)fitness11(x,DATA);%%�ۺ�����ָ��2
% fitness = @(x)fitness12(x,DATA);%% ��߶�������
% fitness = @(x)fitness13(x,DATA);%% ��߶�������
% fitness = @(x)fitness14(x,DATA);%% ��߶�ģ����
%% �Ż�Ѱ��
SearchAgents_no =40;                   % ����
Max_iteration = 10;                    % ����������
dim = 2;                               % �Ż���������
lb = [1000, 8];                 % ����ȡֵ�½�(a,k)
ub = [4000, 15];                 % ����ȡֵ�Ͻ�
%% ����Ѱ��
[Best_score,Best_pos,Convergence_curve,ak]=PSO(SearchAgents_no,Max_iteration,lb ,ub,dim,fitness);
%% ��ʾ���Ų���������VMD
disp(['����KֵΪ��',num2str(round(Best_pos(1,2)))])
disp(['����alphaֵΪ��',num2str(round(Best_pos(1,1)))])
disp(['����ָ��Ϊ��',num2str(Best_score)])
[u, u_hat, omega] = VMD(DATA.data, round(Best_pos(1,1)), DATA.tau, round(Best_pos(1,2)), DATA.dc, DATA.init, DATA.tol);
%% ���ӻ�
%% ��Ӧ������
figure
plot(Convergence_curve,'linewidth',1.5);
title('��������')
xlabel('��������')
ylabel('��Ӧ��ֵ')
xlim([1,Max_iteration]);
grid on
saveas(gcf,['���ͼ/1',],'fig');
[m,n]=size(u);
%% �ֽ������ӻ�
figure
subplot(m+1,1,1);
plot(DATA.data,'k');grid on;
title('ԭʼͼ��');
for i = 1:m
    subplot(m+1,1,i+1);
    plot(u(i,:),'k');
    title(['PSO-VMD�ֽ��ź�',num2str(i),':']);
end
saveas(gcf,['���ͼ/2',],'fig');
set(gcf,'position',[50 50 650 650])%%����������С
%% Ƶ��ͼ���ӻ�
figure
for i=1:m
subplot(m,1,i)
%% FFT �任
[cc,y_f]=hua_fft(u(i,:),fs,1);
plot(y_f,cc,'b','LineWIdth',1.5);
% hua_fft(u(i,:),fs,1)
ylabel(['FFT of IMF',num2str(i)]);
axis tight
end
saveas(gcf,['���ͼ/3',],'fig');
%% ���Ų����仯����
figure
plot(round(ak(:,1)),'k-s ','linewidth',1);
title('�ͷ������Ż���������')
xlabel('��������')
ylabel('�ͷ�����')
xlim([1,Max_iteration]);
ylim([lb(:,1),ub(:,1)]);
saveas(gcf,['���ͼ/4',],'fig');
figure
plot(round(ak(:,2)),'k-d ','linewidth',1);
title('�ֽ�ģ̬���Ż���������')
xlabel('��������')
ylabel('�ֽ�ģ̬��')
xlim([1,Max_iteration]);
ylim([lb(:,2),ub(:,2)]);
saveas(gcf,['���ͼ/5',],'fig');
%% �ֽ���3D��ͼ
u1=u(1:round(Best_pos(1,2)),:);
plot3imf(u1);%%Ϊ����Ժ������ã�������plot3imf�ļ�
title('�ֽ���3D��ͼ')
saveas(gcf,['���ͼ/6',],'fig');

save ���ͼ/res
