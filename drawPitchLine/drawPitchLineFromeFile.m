function drawpitch2LineFromeFile(filename,whichPicture,id)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

max2Time = 10.0;

[times2_ms  pitch2]=tex2tread(filename,'%f\t%f');

%分割文件 s
total2 = length(times2_ms);
duration2 = times2_ms(length(times2_ms)) - times2_ms(1);
x2 = duration2 / (total2 - 1);
count2 = ceil(duration2/max2Time);



for i=1:count2
    if(i==whichPicture)
        startSample2 = floor(((i-1)*max2Time -  times2_ms(1))/x2)+1;
        endSample2 = floor(i*max2Time - times2_ms(1))/x2;
        if startSample2 <= 0 
            startSample2 = 1;
        end
        if endSample2 > length(times2_ms)
            endSample2 = length(times2_ms);
        end
        tarr = times2_ms(startSample2:endSample2);
        parr = pitch2(startSample2:endSample2);
        str = [filename num2str(i)];
%         figure('NumberTitle', 'off', 'Name', str);
        subplot(2,1,id);
        plot(tarr,parr,'*r');
    end
end

end

