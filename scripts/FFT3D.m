function Image = FFT3D(Kspace)

Image =  ifftshift(ifftshift(ifftshift(   fft( fft( fft(   fftshift(fftshift(fftshift(Kspace,1),2),3),[],1),[],2),[],3),3),2),1);
