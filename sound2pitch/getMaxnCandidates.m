function [maxnCandidates] = getMaxnCandidates(arr)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
maxnCandidates = 0;
for i=1:length(arr)
    if(arr(i) > maxnCandidates)
        maxnCandidates = arr(i);
    end
end
end

