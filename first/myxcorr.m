function r=myxcorr(x)
%这个函数假定，窗大小是1536，x长度必须是偶数
r=zeros(1,length(x));
windows=1536;
tmax = windows;
for t=0:length(x)/2-1
for i=1:length(x)/2
    if(t<=tmax)
        r(t+1) = r(t+1) +  x(i)*x(i+t)*(1-t/tmax); 
    else
        r(t+1) = 0;
    end
    
end
end
end