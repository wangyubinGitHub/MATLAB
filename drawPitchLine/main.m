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
filename = 'W.txt';
filename1 = 'W.log';

whichPicture = 1;
figure;
drawPitchLineFromeFile(filename,whichPicture,1);
drawPitchLineFromeFile(filename1,whichPicture,2);
end

