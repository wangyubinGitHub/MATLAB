function [result] = fftcorr(arr)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
   arrlen = length(arr);
   len = 2*arrlen;
   result = ifft(conj(fft(arr,len)).*(fft(arr,len)));
   %result = ifftshift(result);
   result = [result(arrlen+2:length(result)),result(1:arrlen)];%��Ϊfft��0���µ�0��ƫ�Ƶ�����
end

