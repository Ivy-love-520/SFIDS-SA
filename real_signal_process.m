clear;clc;close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
fs = 5120; 
data = load('E01_240730_091352386.txt');
 transposedData = data';
dlmwrite('data.txt', transposedData);
 signal = load('data.txt');
s2 = signal(3,:);
Sig=s2;
Os = s2;
AM = abs(fft(Os));
sample = length(Os);
dt = 0.002;
t = 0.002:dt:dt*sample;
df = 1/(dt*sample);         
f = 0:df:(sample - 1)*df;     
freq = 250;  
[~,id] = min(abs(f - freq));   
TFRgst = TFT(s2,[],dt,'GST',0.2,20);
DATA.data = Sig;
DATA.tau = 0;
DATA.dc = 0;
DATA.init = 1;
DATA.tol = 1e-7;
fitness = @(x)fitness1(x,DATA); 
SearchAgents_no = 10;                 
Max_iteration = 10;                    
dim = 2;                               
lb = [100, 2];                 
ub = [3000, 10];               
[Best_score,Best_pos,Convergence_curve,ak]=PSO(SearchAgents_no,Max_iteration,lb,ub,dim,fitness);
disp(['The best K：',num2str(round(Best_pos(1,2)))])
disp(['The best alpha：',num2str(round(Best_pos(1,1)))])
disp(['The best：',num2str(Best_score)])
K = 7;
[imf,residual] = vmd(Sig,"NumIMFs",K);
windowLength = 256; 
overlap = 0.5; 
nfft = 512; 
for i = 1:K
    eval(['imf', num2str(i), '= imf(:,', num2str(i), ');']);
    [S,F,T] = spectrogram(eval(['imf', num2str(i)]), windowLength, floor(windowLength * overlap), nfft, fs);
    eval(['STFT_imf', num2str(i), ' = S;']);
    eval(['Freq_STFT_imf', num2str(i), ' = F;']);
    eval(['Time_STFT_imf', num2str(i), ' = T;']);
end
for i = 1:K
    eval(['AM_imf', num2str(i), '= abs(fft(', 'imf', num2str(i), '));']);
end
figure(1);
subplot(2,1,1);
plot(t,s2,'b');xlim([0 3]);xlabel('Time(s)');ylabel('Amplitude');
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
    eval(['plot(f(1,1:id),AM_imf', num2str(i), '(1:id,1));']);
    xlabel('Frequency(Hz)');
    ylabel('Amplitude');
end
normalized_AM = AM / max(AM);
area_AM_imf_array = zeros(1, K); 
for i = 1:K
    eval(['normalized_AM_imf', num2str(i), '= AM_imf', num2str(i), ' / max(AM_imf', num2str(i), ');']);
    eval(['area_AM_imf', num2str(i), '= trapz(f, normalized_AM_imf', num2str(i), ');']);
    area_AM_imf_array(i) = eval(['area_AM_imf', num2str(i)]);
end
disp('normal AM area：');
disp(trapz(f, normalized_AM));
for i = 1:K
    disp(['AM_imf', num2str(i), 'normal area：', num2str(area_AM_imf_array(i))]);
end
threshold = 0.1; 
over_decomposition = false; 
for i = 1 : K  
      if area_AM_imf_array(i) < threshold
        over_decomposition = true;
        break;
      end
end
if over_decomposition
    disp('over-decomposition');
else
    disp('under-decomposition');
end
max_values = zeros(1, K);  
counted=0;
for i = 1:K  
    eval(['AM_imf_data = abs(fft(', 'imf', num2str(i), '));']);  
    AM_imf_data = eval(['AM_imf', num2str(i)]);  
    [max_value, max_index] = max(AM_imf_data);  
    max_values(i) = max_value; 
    if max_index == 1  
        ratio = max_value / AM_imf_data(max_index + 1);  
    elseif max_index == length(AM_imf_data)  
        ratio = max_value / AM_imf_data(max_index - 1);  
    else  
        left_ratio = max_value / AM_imf_data(max_index - 1);  
        right_ratio = max_value / AM_imf_data(max_index + 1);  
        ratio = max([left_ratio, right_ratio]);  
    end  
    disp(['IMF', num2str(i), ' The ratio of the maximum value in the frequency domain to the values on both sides：', num2str(ratio)]);  
    area_ratio = ratio / area_AM_imf_array(i);  
    disp(['IMF', num2str(i), ' The ratio of the proportion to the area：', num2str(area_ratio)]);  
    if area_ratio > 0.5  
        disp(['IMF', num2str(i), ' Satisfies the condition that the ratio of the proportion to the area is greater than 0.5']);  
        counted = counted + 1;  
        if max_index <= 5  
            neighbor_avg = mean(AM_imf_data(1:min(max_index + 10, length(AM_imf_data))));  
        elseif max_index >= length(AM_imf_data) - 5  
            neighbor_avg = mean(AM_imf_data(max_index - 10:end));  
        else  
            left_avg = mean(AM_imf_data(max_index - 5:max_index - 1));  
            right_avg = mean(AM_imf_data(max_index + 1:min(max_index + 5, length(AM_imf_data))));  
            neighbor_avg = mean([left_avg, right_avg]);  
        end 
        AM_imf_data(AM_imf_data > neighbor_avg*0.5) = neighbor_avg*0.1; 
        windowSize = 5;  
        AM_imf_data = movmean(AM_imf_data, windowSize);
        eval(['AM_imf', num2str(i), ' = AM_imf_data;']);    
        original_imf = eval(['imf', num2str(i)]);  
        updated_imf = ifft(AM_imf_data, 'symmetric'); 
        eval(['imf', num2str(i), ' = updated_imf;']);
    end  
end  

figure(4); 
for i = 1:K  
    eval(['AM_imf_data = abs(fft(', 'imf', num2str(i), '));']);  
    subplot(K, 1, i);  
    plot((f(1:id)), AM_imf_data(1:id),'b'); 
    xlabel('Frequency (Hz)');  
    ylabel('Amplitude');  
    title(['IMF Frequency Domain - IMF', num2str(i)]);  
end  
sgtitle('Updated IMF Frequency Domain');  
IMF = zeros(length(Os), K);  
figure(6);
for i = 1:K
    current_AM_imf = eval(['AM_imf', num2str(i)]);
    Z=fft(s2);
    Y=(current_AM_imf)'.*exp(1i*angle(Z));
    imf=real(ifft(Y)); 
    IMF(:, i) = imf;  
    subplot(K, 1, i);
    plot(t, imf,'b'); 
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Inverse Fourier Transform of AM_imf', num2str(i)]);
end
imf_sum=0;
for i = 1:K
    imf_sum = IMF(:,i) + imf_sum;
end

%% Perform K-SVD
NoiseSigma = NoiseEstimate(imf_sum);
Params.n = 128;
Params.m = Params.n * 2;
Params.E = 20 * NoiseSigma * getConstant(Params.n);
tic
[y_KSVD] = CleanKSVD(imf_sum, Params);
toc
AM = abs(fft(y_KSVD));
figure;
subplot(2,1,1)
plot(t, y_KSVD, 'b');
subplot(2,1,2);
plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
TFRgst_y = TFT(y_KSVD,[],dt,'GST',0.2,20);
AM_s2 = abs(fft(s2));
figure(11)
subplot(6,2,1);
plot(t,s2,'b');xlabel('Time(s)');ylabel('Amplitude');xlim([0,3]);
temp = abs(TFRgst(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(6,2,3);
plot((f(1:id)),AM_s2(1:id),'b');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
subplot(6,2,5);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');xlim([0,3]);
set(gca,'ydir','normal');
colormap jet
subplot(6,2,2);
plot(t,y_KSVD,'b');xlabel('Time(s)');ylabel('Amplitude');xlim([0,3]);
temp = abs(TFRgst_y(:,1:id)');
tfr = temp ./ max(max(temp));
subplot(6,2,4);
plot((f(1:id)),AM(1:id),'b');xlabel('Frequency(Hz)');ylabel('Amplitude');
subplot(6,2,6);
imagesc(t,f(1:id),tfr);xlabel('Time(s)');ylabel('Frequency(Hz)');xlim([0,3]);
set(gca,'ydir','normal');
colormap jet
