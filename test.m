clear;clc;close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%合成信号
%% 加载txt格式的数据
s_clean = load('s_clean.txt');        %s_clean为原始干净信号
s_clean = (s_clean)';
s1 = load('signal.txt');              %signal为加单频干扰的信号
s2 = load('s_noise.txt');
%% 信号命名
s2 = s1';
% 定义噪声参数
mean_noise = 0;               % 噪声的均值（可以根据需要调整）
std_dev_noise = 0.004;          % 噪声的标准差（可以根据需要调整）
% 生成随机噪声
noise = std_dev_noise * randn(size(s2)) + mean_noise;
% 将噪声添加到 s2
s2 = s2 + noise;                 %s2为加了噪声和单频干扰的信号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig = s2;
Os = s2;
fs = 5120; % 采样频率，用于频谱图
AM = abs(fft(Os));
sample = length(Os);
dt = 0.002; 
t = 0.002:dt:dt*sample;
df = 1/(dt*sample);         % 频率间隔
f = 0:df:(sample - 1)*df;     % 频率采样点数
freq = 250;  % 显示时的 最大频率
[~,id] = min(abs(f - freq));   % 求出最大频率时对应的点数位置
TFRgst = TFT(s2,[],dt,'GST',0.2,20);
DATA.data = Sig;
DATA.tau = 0;
DATA.dc = 0;
DATA.init = 1;
DATA.tol = 1e-7;
%% 更换适应度函数
fitness = @(x)fitness1(x,DATA); % 综合评价指标1
%% 优化寻参
SearchAgents_no = 10;                   % 数量
Max_iteration = 10;                    % 最大迭代次数
dim = 2;                               % 优化参数个数
lb = [100, 7];                 % 参数取值下界(a,k)
ub = [3000, 7];                 % 参数取值上界
%% 参数寻参
tic
[Best_score,Best_pos,Convergence_curve,ak]=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fitness);
toc
%% 显示最优参数
disp(['最优K值为：',num2str(round(Best_pos(1,2)))])
disp(['最优alpha值为：',num2str(round(Best_pos(1,1)))])
disp(['最优指标为：',num2str(Best_score)])
%将最优参数代入VMD
K = round(Best_pos(1,2));
[imf,residual] = vmd(Sig,"NumIMFs",K);
% 设置STFT参数
windowLength = 256; % 窗口长度，可根据信号特点调整
overlap = 0.5; % 窗口重叠比例，取值范围0到1
nfft = 512; % FFT点数，通常选择大于等于窗口长度的2的幂次方
% 对每个IMF进行短时傅里叶变化（STFT）
for i = 1:K
    eval(['imf', num2str(i), '= imf(:,', num2str(i), ');']);
    % 进行STFT
    [S,F,T] = spectrogram(eval(['imf', num2str(i)]), windowLength, floor(windowLength * overlap), nfft, fs);
    % 存储STFT结果
    eval(['STFT_imf', num2str(i), ' = S;']);
    eval(['Freq_STFT_imf', num2str(i), ' = F;']);
    eval(['Time_STFT_imf', num2str(i), ' = T;']);
end
for i = 1:K
    eval(['AM_imf', num2str(i), '= abs(fft(', 'imf', num2str(i), '));']);
end
figure(1);
subplot(2,1,1);
plot(t,s2,'b');xlabel('Time(s)');ylabel('Amplitude');
subplot(2,1,2);
plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
figure(2);
for i = 1:K
    subplot(K, 1, i);
    plot(t, imf(:, i),'b');
    xlabel('Time(s)');
    ylabel('Amplitude');
end
figure(3);
for i = 1:K
    subplot(K, 1, i);
    eval(['plot(f(1,1:id), AM_imf', num2str(i), '(1:id,1), ''b'');']);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 对频率域进行归一化处理
normalized_AM = AM / max(AM);
area_AM_imf_array = zeros(1, K); % 先创建一个空数组来 存储每个area_AM_imf的值
for i = 1:K
    eval(['normalized_AM_imf', num2str(i), '= AM_imf', num2str(i), ' / max(AM_imf', num2str(i), ');']);
    eval(['area_AM_imf', num2str(i), '= trapz(f, normalized_AM_imf', num2str(i), ');']);
    % 将每个area_AM_imf的值存储到数组中
    area_AM_imf_array(i) = eval(['area_AM_imf', num2str(i)]);
end
% 输出每个AM_imf对应的面积结果以及AM的面积结果
disp('AM归一化后面积为：');
disp(trapz(f, normalized_AM));
for i = 1:K
    disp(['AM_imf', num2str(i), '归一化后面积为：', num2str(area_AM_imf_array(i))]);
end
% 判断IMF是否存在过分解
threshold = 0.1; % 设定一个阈值，可根据实际情况调整
over_decomposition = false; % 初始化过分解标志为false
for i = 1 : K  
      if area_AM_imf_array(i) < threshold
        over_decomposition = true;
        break;
      end
end
if over_decomposition
    disp('存在IMF过分解情况');
else
    disp('未发现IMF过分解情况');
end
% 初始化一个数组来存储每个 IMF 的最大值  
max_values = zeros(1, K);  
%计数器
counted=0;
for i = 1:K  
    % 计算当前 IMF 的频率域  
    eval(['AM_imf_data = abs(fft(', 'imf', num2str(i), '));']);  
    AM_imf_data = eval(['AM_imf', num2str(i)]);  
    % 找到频率域中的最大值及其索引  0
    [max_value, max_index] = max(AM_imf_data);  
    max_values(i) = max_value;  % 存储最大值 
    % 计算最大值与邻近值的比值  
    if max_index == 1  
        ratio = max_value / AM_imf_data(max_index + 1);  
    elseif max_index == length(AM_imf_data)  
        ratio = max_value / AM_imf_data(max_index - 1);  
    else  
        left_ratio = max_value / AM_imf_data(max_index - 1);  
        right_ratio = max_value / AM_imf_data(max_index + 1);  
        ratio = max([left_ratio, right_ratio]);  
    end  
    disp(['IMF', num2str(i), ' 频率域最大值与两侧值的比例为：', num2str(ratio)]);  
    % 计算比例与对应 IMF 面积的比值  
    area_ratio = ratio / area_AM_imf_array(i);  
    disp(['IMF', num2str(i), ' 比例与面积的比值为：', num2str(area_ratio)]);  
    % 检查条件并执行替换操作  
    if area_ratio > 0.5  
        disp(['IMF', num2str(i), ' 满足比例与面积的比值大于 0.5 的条件']);  
        counted = counted + 1;
          % 计算邻居值的平均（修改为计算最大值左边5个和右边5个的值的平均）  
        if max_index <= 5  
            neighbor_avg = mean(AM_imf_data(1:min(max_index + 10, length(AM_imf_data))));  
        elseif max_index >= length(AM_imf_data) - 5  
            neighbor_avg = mean(AM_imf_data(max_index - 10:end));  
        else  
            left_avg = mean(AM_imf_data(max_index - 5:max_index - 1));  
            right_avg = mean(AM_imf_data(max_index + 1:min(max_index + 5, length(AM_imf_data))));  
            neighbor_avg = mean([left_avg, right_avg]);  
        end 
        %% 替换频域中大于邻居平均的值  
        AM_imf_data(AM_imf_data > neighbor_avg*0.5) = neighbor_avg*0.1; 
        windowSize = 5;  
        AM_imf_data = movmean(AM_imf_data, windowSize);
        eval(['AM_imf', num2str(i), ' = AM_imf_data;']);  % 更新 IMF 数据  
        % 更新原始 IMF（如果需要保持一致）  
        original_imf = eval(['imf', num2str(i)]);  % 获取原始 IMF  
        updated_imf = ifft(AM_imf_data, 'symmetric'); % 进行逆 FFT  
        eval(['imf', num2str(i), ' = updated_imf;']); % 替换更新后的 IMF  
    end  
end  
% 判断信号是否存在单频干扰
if counted==0
    disp('不存在单频干扰');
    error('程序已被强制终止，因为未发现单频干扰。'); % 强制结束程序 
else
    disp('存在单频干扰');
end
% 绘制更新后的 IMF 分量的频率域图片  
figure(4); % 创建新的图像窗口  
for i = 1:K  
    % 计算并绘制更新后的 IMF 的频率域  
    eval(['AM_imf_data = abs(fft(', 'imf', num2str(i), '));'],'b');  % 计算更新后的频率域  
    % 绘图   
    subplot(K, 1, i);  
    %plot(f, AM_imf_data(1:length(f))); % 绘制频率域图  
    plot((f(1:id)), AM_imf_data(1:id),'b');% 绘制频率域图  
    xlabel('Frequency (Hz)');  
    ylabel('Amplitude');  
    title(['IMF Frequency Domain - IMF', num2str(i)]);  
end  
% 添加全局标题  
sgtitle('Updated IMF Frequency Domain');  
% 新建一个图形窗口用于展示逆变换后的结果（这里可根据实际需求选择合适的展示方式和窗口编号等）
% 假设 K 是IMF的个数，t 是时间向量
IMF = zeros(length(Os), K);  % 初始化IMF数组，存储每个IMF分量
figure(6);
% 循环遍历每个IMF分量进行傅里叶逆变换及相关操作
for i = 1:K
    % 获取当前的AM_imf数据（假设之前的数据存到了对应的AM_imf变量中）
    current_AM_imf = eval(['AM_imf', num2str(i)]);
    % 进行傅里叶逆变
    Z=fft(s2);
    Y=(current_AM_imf)'.*exp(1i*angle(Z));
    imf=real(ifft(Y)); %% 求取逆傅里叶变换
    % 存储IMF分量到数组 IMF 的第i列
    IMF(:, i) = imf;  
    % 绘图展示逆变换后的时域波形（可根据实际需求调整展示方式等细节）
    subplot(K, 1, i);
    plot(t, imf,'b'); % 绘制时域波形，取实部展示，ifft结果可能有虚部，按需决定是否只展示实部
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Inverse Fourier Transform of AM_imf', num2str(i)]);
end
imf_sum=0;
for i = 1:K
    imf_sum = IMF(:,i) + imf_sum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform K-SVD
NoiseSigma = NoiseEstimate(imf_sum);
Params.n = 38;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);
tic
[y_KSVD] = CleanKSVD(imf_sum, Params);
toc
% 如果需要，可以重构降噪后的信号
%denoisedSignal = sum(y_KSVD_all, 2);  % 所有降噪后的IMF分量相加，得到重构信号
AM = abs(fft(y_KSVD));
% 绘制重构的降噪信号
figure;
subplot(2,1,1)
plot(t, y_KSVD, 'b');
subplot(2,1,2);
plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
TFRgst_y = TFT(y_KSVD,[],dt,'GST',0.2,20);
%% 原始信号成图 GST
figure(11)
subplot(2,1,1);
plot(t,s2,'b');xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');
% xlabel('时间/s');
% ylabel('频率/Hz');
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet
%% 处理后的信号成图 GST
figure(12)
subplot(2,1,1);
plot(t,y_KSVD,'b');xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst_y(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');
% xlabel('时间/s');
% ylabel('频率/Hz');
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet
%噪声残差剖面
%将s2与y_KSVD的y坐标统一，以便后续的残差剖面计算
zz_s2=max(s2);
y_KSVD_max = max(y_KSVD);
ccc=zz_s2/y_KSVD_max;
y_new=ccc.*y_KSVD;
error=y_new-s2';     %残差剖面
figure(100);
subplot(3,1,1);
plot(s2);
subplot(3,1,2);
plot(y_KSVD);
subplot(3,1,3);
plot(error);
am_s1=abs(fft(s1));
am_y_new=abs(fft(y_new));
%时间对比图
figure(101);
plot(t,s_clean,'b');hold on;
plot(t,y_KSVD,'r');
%频域对比图
figure(102);
plot(am_s1(1:600)); hold on;
plot(am_y_new(1:600),'r');
figure(103);
plot(t,s_clean);
