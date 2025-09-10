clc;clear;close all
%
addpath(genpath(fileparts(mfilename('fullpath'))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = importdata('08.asc');
% Sig = X(2:end,3);
% fs = 12000;
% t = (0:length(Sig)-1)/fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 合成信号
fs = 1000;           % 采样频率 (Hz)
T = 1656 / fs;       % 信号长度（秒）
N = 1656;            % 信号点数
t = (0:N-1) / fs;    % 时间向量
% 生成雷克子波 (Ricker wavelet)
f0 = 20;             % 雷克子波的中心频率 (Hz)
ricker_wavelet = (1 - 2*(pi*f0*(t - t(round(N/2)))).^2) .* exp(-(pi*f0*(t - t(round(N/2)))).^2);
% 生成单频干扰信号 (例如频率为50 Hz的正弦波)
f_interference = 50;     % 干扰信号频率 (Hz)
interference_signal = 0.5 * sin(2*pi*f_interference*t);  % 干扰信号幅度为0.5
% 叠加信号
Sig= ricker_wavelet + interference_signal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t,Sig,'b');
xlabel('t/s');ylabel('幅值')
title('原始含噪信号');

%% Perform K-SVD
NoiseSigma = NoiseEstimate(Sig);

Params.n = 100;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);

tic
[y_KSVD] = CleanKSVD(Sig, Params);
toc

figure;
plot(t,y_KSVD,'b');
xlabel('t/s');ylabel('幅值')
title('KSVD降噪结果');

window = 50;
Nfrebin = 1024;
SampFreq = 12000;

[Spec1,f] = STFT(Sig(:),SampFreq,Nfrebin,window);
a = abs(Spec1);
figure;
imagesc(t,f,abs(Spec1)); 
axis([0 max(t) 0 SampFreq/2]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');	
b=f';

[Spec2,f] = STFT(y_KSVD(:),SampFreq,Nfrebin,window);
figure;
imagesc(t,f,abs(Spec2));
axis([0 max(t) 0 SampFreq/2]);
set(gcf,'Position',[20 100 320 250]);	 
xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',12);
set(gcf,'Color','w');
%%