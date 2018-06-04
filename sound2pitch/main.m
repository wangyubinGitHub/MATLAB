function main()
%MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
% filename='../test.wav';
%  filename='../sourcefile/source1.wav';
  filename = '../sourcefile/tmpf.wav';
minpitch=75;
timeStep=0;
maxpitch = 500;

maxtime = 10.0;
% 1
 [y, fs] = audioread(filename);
 if(fs ~= 44100)
     return;
 end
%1 end

count  = ceil( length(y) / fs /maxtime);

for i=1:count
    startSample = 1+(i-1)*maxtime*fs;
    endSample = (i)*maxtime*fs;
    if(endSample > length(y))
       endSample = length(y); 
    end
    y1 = y(startSample:endSample);
    sound_to_pitch(y1,timeStep,minpitch,maxpitch);
end
end

