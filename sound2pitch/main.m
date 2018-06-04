function main()
%MAIN 此处显示有关此函数的摘要
%   此处显示详细说明
% filename='../test.wav';
% filename='../sourcefile/source1.wav';
filename = '../sourcefile/tmpf.wav';
minpitch=75;
timeStep=0;
maxpitch = 500;

% 1
 [y, fs] = audioread(filename);
 if(fs ~= 44100)
     return;
 end
%1 end

sound_to_pitch(y,timeStep,minpitch,maxpitch);

end

