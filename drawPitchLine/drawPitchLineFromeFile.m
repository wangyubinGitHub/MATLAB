function drawPitchLineFromeFile(filename)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

[times_ms  pitch]=textread(filename,'%d\t%d');

plot(times_ms,pitch,'*r');

end

