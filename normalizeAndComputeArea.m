% 定义归一化和面积计算函数
function [normalizedData, area] = normalizeAndComputeArea(data, fs)
    % FFT 变换获取频谱
    [cc,y_f]=hua_fft(data,fs,1);
    % 归一化
    normalizedData = cc / max(cc);
    % 计算面积（这里假设频率间隔均匀，直接用积分近似为求和）
    area = sum(normalizedData) * (y_f(2)-y_f(1));
end
