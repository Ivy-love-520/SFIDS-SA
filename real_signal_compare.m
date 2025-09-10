clear;clc;close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %合成信号
%s_clean = load('s_clean.txt');        %s_clean为原始干净信号
% s_clean = (s_clean)';
% s1 = load('signal.txt');              %signal为加单频干扰的信号
% s2 = load('s_noise.txt');
% s2 = s1';
% % 定义噪声参数
% mean_noise = 0;               % 噪声的均值（可以根据需要调整）
% std_dev_noise = 0.004;          % 噪声的标准差（可以根据需要调整）
% % 生成随机噪声
% noise = std_dev_noise * randn(size(s2)) + mean_noise;
% % 将噪声添加到 s2
% s2 = s2 + noise;                 %s2为加了噪声和单频干扰的信号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 导入数据
fs = 5120; % 采样频率，用于频谱图
% 读取纯数值数据
%数据1
data = load('E01_240731_104328334.txt');
%数据2
%data = load('E01_240730_091352386.txt');
% 转置数据
transposedData = data';
% 保存转置后的数据到新文件
dlmwrite('data.txt', transposedData);
signal = load('data.txt');
%设置第几列的信号num
num = 3;
s2 = signal(num,:);
YSJ=signal(num,:);
Sig1=signal(num,:);
Sig2=signal(num,:);
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
K = 8;
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
imf_VMD=imf2+imf3+imf4+imf5+imf6+imf7+imf8;
AM_VMD=abs(fft(imf_VMD));
imfevv=imf1;
for i = 1:K
    eval(['AM_imf', num2str(i), '= abs(fft(', 'imf', num2str(i), '));']);
end
figure;
subplot(2,1,1);
plot(t,s2,'b');xlabel('Time(s)');ylabel('Amplitude');
subplot(2,1,2);
plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
figure;
for i = 1:K
    subplot(K, 1, i);
    plot(t, imf(:, i),'b');
    xlabel('Time(s)');
    xlim([0 3]);
    ylabel('Amplitude');
end
figure;
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
figure; % 创建新的图像窗口  
for i = 1:K  
    % 计算并绘制更新后的 IMF 的频率域  
    eval(['AM_imf_data = abs(fft(', 'imf', num2str(i), '));'],'b');  % 计算更新后的频率域  
    % 绘图   
    subplot(K, 1, i);  
    %plot(f, AM_imf_data(1:length(f))); % 绘制频率域图  
    plot((f(1:id)), AM_imf_data(1:id),'b');% 绘制频率域图  
    xlabel('Frequency (Hz)');  
    xlim([0 150]);
    ylabel('Amplitude');  
    title(['IMF Frequency Domain - IMF', num2str(i)]);  
end  
% 添加全局标题  
sgtitle('Updated IMF Frequency Domain');  
% 新建一个图形窗口用于展示逆变换后的结果（这里可根据实际需求选择合适的展示方式和窗口编号等）
% 假设 K 是IMF的个数，t 是时间向量
IMF = zeros(length(Os), K);  % 初始化IMF数组，存储每个IMF分量
figure;
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
    xlim([0 3]);
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
figure;
plot(t,y_KSVD,'b');
% 如果需要，可以重构降噪后的信号
%denoisedSignal = sum(y_KSVD_all, 2);  % 所有降噪后的IMF分量相加，得到重构信号
AM = abs(fft(y_KSVD));
% % 绘制重构的降噪信号
% figure;
% subplot(2,1,1)
% plot(t, y_KSVD, 'b');xlim([0,3]);xlabel('Time(s)');ylabel('Amplitude');
% subplot(2,1,2);
% plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
TFRgst_y = TFT(y_KSVD,[],dt,'GST',0.2,20);
% TFRgst_s_clean = TFT(s_clean,[],dt,'GST',0.2,20);
%     % 假设 fs 是采样频率，N 是信号长度，t 是时间向量，f 是频率向量
%     fs = 500; % 示例采样频率
%     N = length(TFRgst_s_clean); % 假设 TFRgst_s_clean 是广义 S 变换的结果
%     t = (0:N - 1) / fs; % 时间向量
%     f = (0:N - 1) * (fs / N); % 频率向量
% 
%     % 找到 200Hz 对应的索引
%     max_freq = 150;
%     id = find(f <= max_freq, 1, 'last');
% %% 原始信号成图 GST(时频图)
% figure;
% subplot(2,1,1);
%plot(t,s_clean,'b');xlabel('Time(s)');ylabel('Amplitude');xlim([0,3]);
% temp = abs(TFRgst_s_clean(:,1:id)');
% tfr = temp ./ max(max(temp));
% subplot(2,1,2);
% imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');xlim([0,3]);
% % xlabel('时间/s');
% % ylabel('频率/Hz');
% set(gca,'ydir','normal');
% % set(gca,'YTick',0:50:200);
% colormap jet
%% 加噪信号成图 GST(时频图)
figure;
subplot(2,1,1);
plot(t,s2,'b');xlabel('Time(s)');ylabel('Amplitude');xlim([0,3]);ylim([-0.1,0.1]);
temp = abs(TFRgst(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');xlim([0,3]);
% xlabel('时间/s');
% ylabel('频率/Hz');
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet
%% 处理后的信号成图 GST(时频图)
AM_y_KSVD=abs(fft(y_KSVD));
figure;
subplot(3,1,1);
plot(t,y_KSVD,'b');xlim([0 3]);xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst_y(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(3,1,2);
plot((f(1:id)),AM(1:id),'b');xlim([0 150]);xlabel('Frequency(Hz)');ylabel('Amplitude');
subplot(3,1,3);
imagesc(t,f(1:id),tfr);xlim([0 3]);xlabel('Time(s)');ylabel('Frequency(Hz)');
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
%加噪信号，本文方法去噪以及残差剖面
    figure;
    subplot(3,1,1)
    plot(t,s2,'b');
    ax = gca;
    xlim([0,3]);
    ylim([-0.1,0.1]);
    xlabel('Time(s)');
    ylabel('Amplitude');
    set(ax, 'FontSize', 18);
    subplot(3,1,2)
    plot(t,y_KSVD, 'b');
    ax = gca;
    xlim([0,3]);
    ylim([-0.1,0.1]);
    xlabel('Time(s)');
    ylabel('Amplitude')
    set(ax, 'FontSize', 18);
    subplot(3,1,3);
    plot(t,error,'b');
    ax = gca;
    xlim([0,3]);
    ylim([-0.1,0.1]);
    xlabel('Time(s)');
    ylabel('Amplitude');
    set(ax, 'FontSize', 18);
%am_s1=abs(fft(s1));
am_y_new=abs(fft(y_new));
%时间对比图
% figure;
% plot(t,s_clean,'b');hold on;xlim([0,3]);
% plot(t,y_KSVD,'r');
%频域对比图
% figure;
% plot(am_s1(1:600)); hold on;
% plot(am_y_new(1:600),'r');
% figure;
% plot(t,y_KSVD,'b');
% figure;
% plot(t,s_clean,'b');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 计算信噪比
% signal_power = mean(s_clean.^2);    % 干净信号功率
% % 计算去噪后与干净信号的误差（实际残余噪声）的功率
% noise_y_KSVD = (s_clean)' - y_KSVD; % 计算噪声
% noise_power_y_KSVD = mean(noise_y_KSVD.^2); % 噪声的功率
% % 计算信噪比
% snr_y_KSVD = 10 * log10(signal_power / noise_power_y_KSVD);
% % 输出结果
% fprintf('SNR after denoising (y_KSVD): %.2f dB\n', snr_y_KSVD);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 计算信号的结构相似性（SSIM）
% [ssim_value, ssim_map] = ssim((y_KSVD)', s_clean);
% fprintf('SSIM between y_KSVD and s_clean: %.4f\n', ssim_value);    % 输出结果
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 计算归一化互相关 (NCC)
% y_KSVD = (y_KSVD)';
% numerator = sum((s_clean - mean(s_clean)) .* (y_KSVD - mean(y_KSVD))); % 交叉协方差
% denominator = sqrt(sum((s_clean - mean(s_clean)).^2) * sum((y_KSVD - mean(y_KSVD)).^2)); % 各自方差乘积的平方根
% ncc_value = numerator / denominator; % 归一化互相关
% fprintf('NCC after KSVD denoising: %.4f\n', ncc_value);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小波变换等方法
%% 载入噪声信号数据，并且和程序放置在同一个文件夹下
% signal = load('s_noise.txt');
% s_clean = load('s_clean.txt');
% s_clean = (s_clean)';
% YSJ= signal;
t1=clock;
fs = 5120; % 采样频率，用于频谱图
sample = length(YSJ);
dt = 0.002;
t = 0.002:dt:dt*sample;
df = 1/(dt*sample);           % 频率间隔
f = 0:df:(sample - 1)*df;     % 频率采样点数
freq = 250;                   % 显示时的 最大频率
[~,id] = min(abs(f - freq));  % 求出最大频率时对应的点数位置
dt = 0.002;
addpath(genpath(fileparts(mfilename('fullpath'))));
[c,l]=size(YSJ);
Y=[];
for i=1:c
    Y=[Y,YSJ(i,:)];
end
[c1,l1]=size(Y);
X=[1:l1];
%% 绘制噪声信号图像
figure;
plot(X,Y,'b');
xlabel('Time(s)');
ylabel('Amplitude');
title('原始信号');
%% 硬阈值处理
lev=2;
xd=wden(Y,'heursure','h','one',lev,'db10');%硬阈值去噪处理后的信号序列
figure;
plot(X,xd,'b')
xlabel('Time(s)');
ylabel('Amplitude');
title('硬阈值去噪处理')
set(gcf,'Color',[1 1 1])
%% 软阈值处理
lev=2;
xs=wden(Y,'heursure','s','one',lev,'db10');%软阈值去噪处理后的信号序列
figure;
plot(X,xs,'b')
xlabel('Time(s)');
ylabel('Amplitude');
title('软阈值去噪处理')
set(gcf,'Color',[1 1 1])
%% 固定阈值后的去噪处理
lev=10;
xz=wden(Y,'sqtwolog','s','sln',lev,'db10');%固定阈值去噪处理后的信号序列
figure;
plot(X,xz,'b');
xlabel('Time(s)');
ylabel('Amplitude');
title('固定阈值后的去噪处理')
set(gcf,'Color',[1 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2=clock;
tim=etime(t2,t1);
disp(['------------------运行耗时',num2str(tim),'秒-------------------'])
TFRgst_xz = TFT(xz,[],dt,'GST',0.2,20);
%% 处理后的信号成图 GST(时频谱图)
figure;
subplot(2,1,1);
plot(xz,'b');xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst_xz(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');
% xlabel('Time(s)');
% ylabel('Frequency(Hz)');
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet
sgtitle('固定阈值的时频谱图');  
TFRgst_xd = TFT(xd,[],dt,'GST',0.2,20);
%% 处理后的信号成图 GST(时频谱图)
figure;
subplot(2,1,1);
plot(t,xd,'b');xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst_xd(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');
% xlabel('Time(s)');
% ylabel('频率/Hz');
sgtitle('硬阈值的时频谱图');  
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%字典学习
addpath(genpath(fileparts(mfilename('fullpath'))));
%X = importdata('08.asc');
% data = load('E01_240730_091352386.txt');
%data = load('E01_240731_104328334.txt');
% s_clean = load('s_clean.txt');
% X = load('s_noise.txt');
Sig1 = X;
fs = 12000;
t = (0:length(Sig1)-1)/fs;
Os = X;
AM = abs(fft(Os));
sample = length(Os);
dt = 0.002;
t = 0.002:dt:dt*sample;
df = 1/(dt*sample);         % 频率间隔
f = 0:df:(sample - 1)*df;     % 频率采样点数
freq = 250;  % 显示时的 最大频率
[~,id] = min(abs(f - freq));   % 求出最大频率时对应的点数位置
figure;
plot(t,Sig1,'b');
xlabel('Time(s)');ylabel('Amplitude')
title('原始含噪信号');
%% Perform K-SVD
NoiseSigma = NoiseEstimate(Sig2);
Params.n = 100;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);
tic
[y_KSVD1] = CleanKSVD(Sig2, Params);
toc
figure;
plot(t,y_KSVD1,'b');
xlabel('Time(s)');ylabel('Amplitude')
title('KSVD降噪结果');


% figure;
% plot(t,s_clean,'b');hold on;
% plot(t,y_KSVD1,'r');

TFRgst_y1 = TFT(y_KSVD1,[],dt,'GST',0.2,20);
%% 处理后的信号成图 GST
figure;
subplot(2,1,1);
plot(t,y_KSVD,'b');xlabel('Time(s)');ylabel('Amplitude');
temp = abs(TFRgst_y1(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(2,1,2);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');
% xlabel('Time(s)');
% ylabel('Frequency(Hz)');
set(gca,'ydir','normal');
% set(gca,'YTick',0:50:200);
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 假设我们有 12 组数据对，这里简单生成一些示例数据
numPairs = 10;
% t = linspace(0, 10, 100); % 定义时间向量


% TFRgst_s_clean = TFT(s_clean,[],dt,'GST',0.2,20);
TFRgst_s2 = TFT(s2,[],dt,'GST',0.2,20);

% % 创建一个新的图形窗口
% figure;
%     subplot(5, 2, 1);
%     % 绘制曲线
%     plot(t,s_clean, 'b');
%     title('合成信号');
%     xlabel('Time(s)');
%     ylabel('Amplitude');
%     xlim([0 3]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 假设 fs 是采样频率，N 是信号长度，t 是时间向量，f 是频率向量
%     fs = 500; % 示例采样频率
%     N = length(TFRgst_s_clean); % 假设 TFRgst_s_clean 是广义 S 变换的结果
%     t = (0:N - 1) / fs; % 时间向量
%     f = (0:N - 1) * (fs / N); % 频率向量
%     % 找到 200Hz 对应的索引
%     max_freq = 150;
%     id = find(f <= max_freq, 1, 'last');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(5, 2, 2);
%     temp = abs(TFRgst_s_clean(:,1:id)');
%     tfr = temp ./ max(max(temp));
%     imagesc(t,f(1:id),tfr);
%     %xlabel('Time(s)');ylabel('Frequency(Hz)');
%     xlim([0 3]);
%     xlabel('Time(s)');
%     ylabel('Frequency(Hz)');
%     set(gca,'ydir','normal');
%     % set(gca,'YTick',0:50:200);
%     colormap jet
    figure;
    subplot(4, 2, 1);
    % 绘制曲线
    plot(t,s2, 'b');xlim([0,3])
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontSize', 18);
    subplot(4, 2, 2);
    temp = abs(TFRgst_s2(:,1:id)');
    tfr = temp ./ max(max(temp));
    imagesc(t,f(1:id),tfr);
    %xlabel('Time(s)');ylabel('Frequency(Hz)');
    xlim([0 3]);
    set(gca,'ydir','normal');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Frequency(Hz)', 'FontName', 'Times New Roman', 'FontSize', 18);
    % set(gca,'YTick',0:50:200);
    colormap jet
    TFRgst_imf_VMD = TFT(imf_VMD,[],dt,'GST',0.2,20);
    subplot(4, 2, 3);
    % 绘制曲线
    plot(t,xz, 'b');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontSize', 18);
    xlim([0 3]);
    subplot(4, 2, 4);
    temp = abs(TFRgst_xz(:,1:id)');
    tfr = temp ./ max(max(temp));
    imagesc(t,f(1:id),tfr);
    %xlabel('Time(s)');ylabel('Frequency(Hz)');
    xlim([0 3]);
    set(gca,'ydir','normal');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Frequency(Hz)', 'FontName', 'Times New Roman', 'FontSize', 18);
    % set(gca,'YTick',0:50:200);
    colormap jet
    subplot(4, 2, 5);
    % 绘制曲线
    plot(t,y_KSVD1, 'b');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontSize', 18);
    xlim([0 3]);
    subplot(4, 2, 6);
    temp = abs(TFRgst_y1(:,1:id)');
    tfr = temp ./ max(max(temp));
    imagesc(t,f(1:id),tfr);
    %xlabel('Time(s)');ylabel('Frequency(Hz)');
    xlim([0 3]);
    set(gca,'ydir','normal');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Frequency(Hz)', 'FontName', 'Times New Roman', 'FontSize', 18);
    % set(gca,'YTick',0:50:200);
    colormap jet
    subplot(4,2, 7);
    % 绘制曲线
    plot(t,y_KSVD, 'b');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontSize', 18);
    xlim([0 3]);
    subplot(4, 2, 8);
    temp = abs(TFRgst_y(:,1:id)');
    tfr = temp ./ max(max(temp));
    imagesc(t,f(1:id),tfr);
    %xlabel('Time(s)');ylabel('Frequency(Hz)');
    xlim([0 3]);
    set(gca,'ydir','normal');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    xlabel('Time(s)', 'FontName', 'Times New Roman', 'FontSize', 18);
    ylabel('Frequency(Hz)', 'FontName', 'Times New Roman', 'FontSize', 18);
    % set(gca,'YTick',0:50:200);
    colormap jet