function Image = FFT2D(Kspace)
Image =  ifftshift(ifftshift(   fft( fft(   ifftshift(ifftshift(Kspace,1),2),[],1),[],2),2),1);
