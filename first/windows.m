function [ r ] = windows( arr )
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
win=hann(length(arr));
r=arr.*win;
end

