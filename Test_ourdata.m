clc;clear;close all
%
addpath(genpath(fileparts(mfilename('fullpath'))));
%%
X = importdata('08.asc');

Sig = X(2:end,3);
fs = 12000;
t = (0:length(Sig)-1)/fs;

figure;
plot(t,Sig,'b');
xlabel('t/s');ylabel('幅值')
title('原始含噪信号');

%% Perform K-SVD
NoiseSigma = NoiseEstimate(Sig);

Params.n = 200;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);

tic
[y_KSVD] = CleanKSVD(Sig, Params);
toc

figure;
plot(t,y_KSVD,'b');
xlabel('t/s');ylabel('幅值')
title('KSVD降噪结果');

window = 128;
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