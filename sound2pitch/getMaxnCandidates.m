function [maxnCandidates] = getMaxnCandidates(arr)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
maxnCandidates = 0;
for i=1:length(arr)
    if(arr(i) > maxnCandidates)
        maxnCandidates = arr(i);
    end
end
end

