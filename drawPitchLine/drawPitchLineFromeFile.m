function drawPitchLineFromeFile(filename)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

[times_ms  pitch]=textread(filename,'%d\t%d');

plot(times_ms,pitch,'*r');

end

