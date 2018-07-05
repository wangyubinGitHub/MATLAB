function main()

x1 = [1,2,3,7,9,8,3,7]';
x2 = [1,2,3,7,9,8,3,7]';

size = length(x1)*2-1;

fft(x1,size)

xcorr(x1,x2)
ifft(fft(x1,size).*conj(fft(x2,size)))
fftshift(ifft(fft(x1,size).*conj(fft(x2,size))))

end

