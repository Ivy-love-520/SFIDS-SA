function [TFR,rec,s_l,A] = TFT(x,r,dt,method,p,lamda)

x = x(:);
n = length(x);
h = fix(n/2) + 1;
sym = 0;
if rem(n,2) == 0
    sym = 1;
end
t = 0:dt:(n-1)*dt;
x_hat = fft(x);
%% windowing
if strcmp(method,'ST')
    s_l = (1:h);
end

if strcmp(method,'GST')
    s_l = (1:h);
    s_l = (lamda*s_l.^p);
end

A = zeros(h);
for  ii = 1:h
    w = g(h,s_l(ii));
    cc = w(h+2-ii:end);
    cc(h+1:end) = [];
    A(ii,:) = cc;
end

A(n,n) = 0;
A(h+1:end,h+1:end) = rot90(A(2:h-sym,2:h-sym),2);
if strcmp(method,'ST')
    A = A';
end
A = A ./ repmat(sum(A,2),1,n); % 确少这一步无法进行准确的信号重构的振幅，波形是一样的，但振幅不一样
%% Computing TFR
TFR = ifft(A.*repmat(x_hat,1,n));
rec = real(sum(TFR,2));
%% 
% figure()
% subplot(211)
% plot(t,rec,'b.',t,x,'k');
% legend('REC','Inp');
% axis tight;
% 
% subplot(212)
% imagesc(t,linspace(0,1/2/dt,h),abs(TFR(:,1:h).'));
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');

end

function gauss = g(m,s)
v = (-m:m-1).^2;
gauss = exp(v*(-1*2*pi^2/s^2));
end

