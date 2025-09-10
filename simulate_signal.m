clc
clear all
% 定义时间采样点数
n = 1000;  
% 时间序列，这里假设采样间隔为0.001，你可以根据实际需求调整
dt = 0.001;  
t = (0 : n - 1) * dt;  
% 雷克子波的主频，单位Hz，可按需调整
fc = 20;  
% 生成雷克子波
w = (1 - 2 * (pi * fc * (t - (n - 1) * dt / 2)).^2).* exp(-(pi * fc * (t - (n - 1) * dt / 2)).^2);
AM_w=abs(fft(w));
% 设置采样率和时间向量  
fs = 1000; % 采样率（Hz）  
t = 0:1/fs:1-1/fs; % 时间向量（持续1秒）  
% 定义信号参数  
high_frequency = 140; % 特别高的频率（Hz）  
low_frequencies = [145.1,145.2,145.3,145.4,145.5,145.6,145.7,145.7,145.8,145.9,146.1,146.2,146.3]; % 相对较低的频率（Hz）  
amplitude_high = 0.1;   % 高频信号的振幅  
% 生成高频正弦波信号  
high_freq_signal = amplitude_high * sin(2 * pi * high_frequency * t);  
% 把雷克子波和单频干扰信号相加，得到最终混合信号
combined_signal = w + high_freq_signal;
s2=combined_signal;



Os = s2;
AM = abs(fft(Os));

sample = length(Os);

dt = 0.002;
t = 0.002:dt:dt*sample;

df = 1/(dt*sample);         % 频率间隔
f = 0:df:(sample - 1)*df;     % 频率采样点数
freq = 250;  % 显示时的 最大频率
[~,id] = min(abs(f - freq));   % 求出最大频率时对应的点数位置
DATA.data = s2;

DATA.tau = 0;
DATA.dc = 0;
DATA.init = 1;
DATA.tol = 1e-7;

%% 更换适应度函数
fitness = @(x)fitness1(x,DATA); % 综合评价指标1
% fitness = @(x)fitness2(x,DATA); % 包络熵
% fitness = @(x)fitness3(x,DATA); % 包络谱峭度
% fitness = @(x)fitness4(x,DATA); % 幅值谱熵
% fitness = @(x)fitness5(x,DATA); % 模糊熵
% fitness = @(x)fitness6(x,DATA); % 皮尔逊相关系数
% fitness = @(x)fitness7(x,DATA); % 峭度值
% fitness = @(x)fitness8(x,DATA); % 样本熵
% fitness = @(x)fitness9(x,DATA); % 排列熵
% fitness = @(x)fitness10(x,DATA); % 信息熵

% fitness = @(x)fitness11(x,DATA); % 综合评价指标2
% fitness = @(x)fitness12(x,DATA); % 多尺度排列熵
% fitness = @(x)ftness13(x,DATA); % 多尺度样本熵
% fitness = @(x)fitness14(x,DATA); % 多尺度模糊熵

%% 优化寻参
SearchAgents_no = 10;                   % 数量
Max_iteration = 10;                    % 最大迭代次数
dim = 2;                               % 优化参数个数
lb = [100, 2];                 % 参数取值下界(a,k)
ub = [3000, 10];                 % 参数取值上界

%% 参数寻参
[Best_score,Best_pos,Convergence_curve,ak]=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fitness);

%% 显示最优参数
disp(['最优K值为：',num2str(round(Best_pos(1,2)))])
disp(['最优alpha值为：',num2str(round(Best_pos(1,1)))])
disp(['最优指标为：',num2str(Best_score)])

%将最优参数代入VMD
K = round(Best_pos(1,2));
[imf,residual] = vmd(s2,"NumIMFs",K);

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
plot(t,s2);xlabel('Time(s)');ylabel('Amplitude');
subplot(2,1,2);
plot((f(1:id)),AM(1:id));xlabel('Frequency(Hz)');ylabel('Amplitude');

figure(5);
for i = 1:K
    subplot(K, 1, i);
    plot(t, imf(:, i));
    xlabel('Time(s)');
    ylabel('Amplitude');
end

figure(6);
for i = 1:K
    subplot(K, 1, i);
    eval(['plot(f(1,1:id),AM_imf', num2str(i), '(1:id,1));']);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude');
end

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
% no_single_frequency_interference = true; % 先假设不存在单频干扰
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
        AM_imf = 0;
    end
end


% %信号重构
 AM_imf_sum = 0; % 根据你的变量类型初始化合适的零值，这里以数值为例，若为矩阵等需要改成zeros初始化相应大小的零矩阵等
% for i = 1:K
%     AM_imf_sum = AM_imf_sum + eval(sprintf('AM_imf%d', i));
% end

AM_imf_sum=AM_imf1+AM_imf2+AM_imf5+AM_imf6+AM_imf7;

% 绘制重构信号的频率域图  
figure(9);  
plot((f(1:id)), AM_imf_sum(1:id));  
xlabel('Frequency(Hz)');  
ylabel('Amplitude');  
title('Frequency Domain of Reconstructed Signal'); 

ZZ=fft(s2);
YY=(AM_imf_sum)'.*exp(1i*angle(ZZ)); 
figure(20);   
plot(AM_imf_sum); 
y=real(ifft(YY)); %% 求取逆傅里叶变换

figure(10); % 创建新的图像窗口用于绘制噪声残差剖面
plot(t,residual); % 绘制噪声残差随时间的变化曲线

% 计算噪声残差
noise_residual = s2 - y;
% 绘制噪声残差剖面
figure(10); % 创建新的图像窗口用于绘制噪声残差剖面
plot(t,noise_residual); % 绘制噪声残差随时间的变化曲线
xlabel('Time(s)');
ylabel('Amplitude');
title('Noise Residual Profile');
% 计算原始信号功率
N = length(s2);
signal_power = sum(s2.^2)/N;
% 计算噪声功率
noise_power = sum(noise_residual.^2)/N;
% 计算信噪比
SNR = 10*log10(signal_power/noise_power);
disp(['信噪比（SNR）为：', num2str(SNR), 'dB']);

% 绘制原始信号与重构信号的对比图
figure(11);
plot(t,s2, 'b',t, y, 'r');
xlabel('Time(s)');
ylabel('Amplitude');
title('Comparison between Original Signal and Reconstructed Signal');
legend('Original Signal', 'Reconstructed Signal');

% 绘制重构信号及重构信号的频率图
figure(12);
subplot(2,1,1);
plot(t,y);
xlabel('Time(s)');
ylabel('Amplitude');
subplot(2,1,2);
plot((f(1:id)), AM_imf_sum(1:id)); 
xlabel('Frequency(Hz)');
ylabel('Amplitude');
