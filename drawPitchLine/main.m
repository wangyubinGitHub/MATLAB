function main()
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

%Ů������
% filename = 'YIN_w.txt';
% filename1 = 'PRT_w.log';

% ��������
% filename = 'YIN_f.txt';
% filename1 = 'PRT_f.log';

% ��������
% filename = 'Z.txt';
% filename1 = 'Z.log';

%�ҵ�����
filename = 'W.txt';
filename1 = 'W.log';

whichPicture = 1;
figure;
drawPitchLineFromeFile(filename,whichPicture,1);
drawPitchLineFromeFile(filename1,whichPicture,2);
end

