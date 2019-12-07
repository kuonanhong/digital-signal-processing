function y = winconv(x,wintype,winamp,winlen)
% Discrere time convolution of a sequence with a window.
% Convolves the sequence X with the
% window specified by wintype having amplitude winamp and length winlen.

% Generate the window
w1 = (window(str2func(wintype),winlen)).'; 
winamp = winamp(:).';
w = winamp.*w1;

% Perform the convolution using FFT
NFFT = 2^(nextpow2(length(x)+winlen));
X = fft(x,NFFT); 
W = fft(w,NFFT);
Y = X.*W;
y = ifft(Y,NFFT);
    