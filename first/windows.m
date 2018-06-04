function [ r ] = windows( arr )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
win=hann(length(arr));
r=arr.*win;
end

