function plot3imf(u)
%   PLOT3IMF 
%   u-----k*N array  k is the number of imf,N is the legnth of the signal

[X,Y] = meshgrid(1:size(u,1),1:size(u,2)); 
Z = u; 
figure;
plot3(Y,X,Z');grid on;

end

