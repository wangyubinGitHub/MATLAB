function r=myxcorr(x)
%��������ٶ�������С��1536��x���ȱ�����ż��
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