clear;clc;close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%% 固定阈值后的去噪处理(固定阈值去噪设置分解层数)
lev=6;
%% 字典个数
n = 128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%合成信号
s_clean = load('s_clean.txt');        %s_clean为原始干净信号
s_clean = (s_clean)';
s1 = load('signal.txt');              %signal为加单频干扰的信号
% s2 = load('s_noise.txt');
% data = xlsread('1.xlsx');
% data = data(:,3);
s2 = s1';
% 定义噪声参数
mean_noise = 0;               % 噪声的均值（可以根据需要调整）
std_dev_noise = 0.005;          % 噪声的标准差（可以根据需要调整）
% 生成随机噪声
noise = std_dev_noise * randn(size(s2)) + mean_noise;
%将噪声添加到 s2
s2 = s2 + noise;                 %s2为加了噪声和单频干扰的信号
Y = s2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s2 = data;
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
%% 固定阈值后的去噪处理
xz=wden(Y,'sqtwolog','s','sln',lev,'db10');%固定阈值去噪处理后的信号序列
figure;
plot(t,xz,'b');
xlabel('Time(s)');
ylabel('Amplitude');
title('固定阈值后的去噪处理')
set(gcf,'Color',[1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算小波固定阈值方法的信噪比
signal_power = mean(s_clean.^2);    % 干净信号功率
% 计算去噪后与干净信号的误差（实际残余噪声）的功率
noise_xz = s_clean- xz; % 计算噪声
noise_power_xz = mean(noise_xz.^2); % 噪声的功率
% 计算信噪比
snr_xz = 10 * log10(signal_power / noise_power_xz);
% 输出结果
fprintf('小波固定阈值SNR: %.4f dB\n', snr_xz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算信号的结构相似性（SSIM）
[ssim_value, ssim_map] = ssim(xz, s_clean);
fprintf('小波固定阈值SSIM: %.4f\n', ssim_value);    % 输出结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算归一化互相关 (NCC)
numerator = sum((s_clean - mean(s_clean)) .* (xz - mean(xz))); % 交叉协方差
denominator = sqrt(sum((s_clean - mean(s_clean)).^2) * sum((xz - mean(xz)).^2)); % 各自方差乘积的平方根
ncc_value = numerator / denominator; % 归一化互相关
fprintf('小波固定阈值NCC: %.4f\n', ncc_value);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform K-SVD
NoiseSigma = NoiseEstimate(Sig);
Params.n = n;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);
tic
[y_KSVD1] = CleanKSVD(Sig, Params);
toc
figure;
plot(t,y_KSVD1,'b');
xlabel('Time(s)');ylabel('Amplitude')
title('KSVD降噪结果');


signal_power = mean(s_clean.^2);    % 干净信号功率
% 计算去噪后与干净信号的误差（实际残余噪声）的功率
noise_y_KSVD1 = (s_clean)'- y_KSVD1; % 计算噪声
noise_power_y_KSVD1= mean(noise_y_KSVD1.^2); % 噪声的功率
% 计算信噪比
snr_y_KSVD1 = 10 * log10(signal_power / noise_power_y_KSVD1);
% 输出结果
fprintf('字典学习SNR: %.4f dB\n', snr_y_KSVD1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算信号的结构相似性（SSIM）
[ssim_value, ssim_map] = ssim((y_KSVD1)', s_clean);
fprintf('字典学习SSIM: %.4f\n', ssim_value);    % 输出结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算归一化互相关 (NCC)
y_KSVD1 = (y_KSVD1)';
numerator = sum((s_clean - mean(s_clean)) .* (y_KSVD1 - mean(y_KSVD1))); % 交叉协方差
denominator = sqrt(sum((s_clean - mean(s_clean)).^2) * sum((y_KSVD1 - mean(y_KSVD1)).^2)); % 各自方差乘积的平方根
ncc_value = numerator / denominator; % 归一化互相关
fprintf('字典学习NCC: %.4f\n', ncc_value);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算加噪信号的信噪比
signal_power = mean(s_clean.^2);    % 干净信号功率
% 计算去噪后与干净信号的误差（实际残余噪声）的功率
noise_s2 = s2 - s_clean; % 计算噪声
noise_power_s2 = mean(noise_s2 .^2); % 噪声的功率
% 计算信噪比
snr_s2  = 10 * log10(signal_power / noise_power_s2 );
% 输出结果
fprintf('含噪信号SNR: %.2f dB\n', snr_s2 );
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算信号的结构相似性（SSIM）
[ssim_value, ssim_map] = ssim(s2, s_clean);
fprintf('含噪信号SSIM: %.4f\n', ssim_value);    % 输出结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算归一化互相关 (NCC)
numerator = sum((s_clean - mean(s_clean)) .* (s2 - mean(s2'))); % 交叉协方差
denominator = sqrt(sum((s_clean - mean(s_clean)).^2) * sum((s2' - mean(s2')).^2)); % 各自方差乘积的平方根
ncc_value = numerator / denominator; % 归一化互相关
fprintf('含噪信号NCC: %.4f\n', ncc_value);
