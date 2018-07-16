function main(start_time,end_time)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

 filename = '/Users/wangyubin/Media/a.txt';
% filename = '/Users/wangyubin/Media/dat.txt';
fid = fopen(filename,'r');

 A = textscan(fid,'%d\t%f');
% A = textscan(fid,'%f\t%f');

fclose(fid);
plot(A{1},A{2},'*r');
return ;
x=A{1};
y=A{2};

startFlag = 1;
endFlag = length(x);
leng = length(x);

for i=1:leng
    if(x(i) >= start_time)
        startFlag = i;
        break;
    end
end

for j=1:leng
    
    if( x(leng - j + 1) <= end_time)
        endFlag = leng - j + 1;
        break;
    end
end
leng
start_time
end_time

startFlag
endFlag

if(startFlag >= endFlag)
   return ; 
end

fx = x(startFlag:endFlag);
fy = y(startFlag:endFlag);
figure;
plot(fx,fy,'*r');


end

