function [result] = fftcorr(arr)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   arrlen = length(arr);
   len = 2*arrlen;
   result = ifft(conj(fft(arr,len)).*(fft(arr,len)));
   %result = ifftshift(result);
   result = [result(arrlen+2:length(result)),result(1:arrlen)];%因为fft补0导致的0点偏移的修正
end

