function main()
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

%女生数据
% filename = 'YIN_w.txt';
% filename1 = 'PRT_w.log';

% 男生数据
% filename = 'YIN_f.txt';
% filename1 = 'PRT_f.log';

% 噪音数据
% filename = 'Z.txt';
% filename1 = 'Z.log';

%我的数据
% filename = 'W.txt';
% filename1 = 'W.log';
%ddd3599f-245b-4fe6-8e70-be33c9d6e4f8
filename = '/Users/wangyubin/Media/316be44c-8fb5-44dd-a0a0-73197c138889.base_pitch';
filename1 = '/Users/wangyubin/Media/316be44c-8fb5-44dd-a0a0-73197c138889.voice_pitch';
whichPicture = 5;
figure;
% drawpitchLineFromeFile(filename,whichPicture,1);
% drawpitchLineFromeFile(filename1,whichPicture,2);

maxTime = 10.0;

[times1_ms  pitch1]=textread(filename,'%f\t%f');

%分割文件 s
total = length(times1_ms);
duration = times1_ms(length(times1_ms)) - times1_ms(1);
x = duration / (total - 1);
count = ceil(duration/maxTime);

[times2_ms  pitch2]=textread(filename1,'%f\t%f');

%分割文件 s
total2 = length(times2_ms);
duration2 = times2_ms(length(times2_ms)) - times2_ms(1);
x2 = duration2 / (total2 - 1);
count2 = ceil(duration2/maxTime);

for i=1:count
    if(i==whichPicture)
        startSample = floor(((i-1)*maxTime -  times1_ms(1))/x)+1;
        endSample = floor((i*maxTime - times1_ms(1))/x);
        if startSample <= 0 
            startSample = 1;
        end
        if endSample > length(times1_ms)
            endSample = length(times1_ms);
        end
        
        startSample2 = floor(((i-1)*maxTime -  times2_ms(1))/x2)+1;
        endSample2 = floor((i*maxTime - times2_ms(1))/x2);
        if startSample2 <= 0 
            startSample2 = 1;
        end
        if endSample2 > length(times2_ms)
            endSample2 = length(times2_ms);
        end
        
        tarr = times1_ms(startSample:endSample);
        parr = pitch1(startSample:endSample);
        str = [filename num2str(i)];
%         figure('NumberTitle', 'off', 'Name', str);

        tarr2 = times2_ms(startSample2:endSample2);
        parr2 = pitch2(startSample2:endSample2);
        str = [filename num2str(i)];
%         figure('NumberTitle', 'off', 'Name', str);

        subplot(4,1,1);
        plot(tarr,parr,'*r');
        
        subplot(4,1,2);
        plot(tarr2,parr2,'*r');
        
        length(tarr2)
        length(tarr)
        subplot(4,1,3);
        plot(tarr,parr,'*r');
        hold on
        plot(tarr2,parr2,'*b');
        
        me = mean(parr2) - mean(parr)
        parr2 = parr2 - me;
                subplot(4,1,4);
        plot(tarr,parr,'*r');
        hold on
        plot(tarr2,parr2,'*b');
    end
end




end

