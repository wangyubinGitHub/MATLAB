function main()
%MAIN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% filename='../test.wav';
%  filename='../sourcefile/source1.wav';
   filename = '../sourcefile/tmpf.wav';
%   filename = '../sourcefile/tmp1.wav';
% filename = '../sourcefile/test1.wav';

minpitch=75;
timeStep=0;
maxpitch = 700;

% maxtime = 10.0;
maxtime = 1000.0
% 1
 [y, fs] = audioread(filename);
 if(fs ~= 44100)
     return;
 end
%1 end

count  = ceil( length(y) / fs /maxtime);
datafile= fopen('data.log', 'a');

for i=1:count
    startSample = 1+(i-1)*maxtime*fs;
    endSample = (i)*maxtime*fs;
    if(endSample > length(y))
       endSample = length(y); 
    end
    y1 = y(startSample:endSample);
%     plot(y1,'r');
    sound_to_pitch(y1,timeStep,minpitch,maxpitch,i,fs*maxtime,datafile);
end

fclose(datafile);
end

