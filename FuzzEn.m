
function   entropy = FuzzEn(x,zr)
% zr=zr*std(x);
[xr,xc] = size(x);
if xr ~= 1	
    x = x';
end
len=length(x);

n=2;
m=2;
N=len-m+1;
for i=1:N
    v(i,:)=x(i:i+m-1);
	mvx(i)=mean(v(i,:));
	vx(i,:)=v(i,:)-mvx(i);
end

B=0;
N1=N-1;  
for i=1:N1
    clear vv;clear maxvv;clear mx1;
	for j=1:m
		vv(1:N1,j)=vx(i,j);
	end
    vv=abs(vv-vx(1:N1,:));
	maxvv=max(vv');
    mx1=exp(-log(2)*(maxvv/zr).^n);   % 指数函数    若为被除数的n次幂型，则改为 
%    mx1=exp(-maxvv.^n/zr);  
    mx1(i)=0;
    cix(i)=sum(mx1);
end
B=sum(cix);

m=m+1;
post= v(2:length(v),m-1);pre= v(1:(length(v)-1),:);vy=[pre post];
for i=1:N1
   mvy(i)=mean(vy(i,:));
   vy(i,:)=vy(i,:)-mvy(i);
end

A=0;
for i=1:N1
    clear vv;clear maxvv;clear my1;
	for j=1:m
		vv(1:N1,j)=vy(i,j);
	end
	vv=abs(vv-vy(1:N1,:));
	maxvv=max(vv');
   my1=exp(-log(2)*(maxvv/zr).^n);   % 指数函数    若为被除数的n次幂型，则改为
% %     my1=exp(-maxvv.^n/zr);  
    my1(i)=0;
    ciy(i)=(sum(my1));
end
A=sum(ciy);

entropy=log(B/A);