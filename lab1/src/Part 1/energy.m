function En = energy(x,wintype,winamp,winlen)
% Short-time energy computation:
% Computes the short-time enery of the sequence x. 
% wintype defines the window type.
% winamp sets the amplitude of the window and 
% winlen is the length of the window in samples

% Enery calculation is achieved by lowpass filtering with window
x2 = x.^2;
En = winconv(x2,wintype,winamp,winlen);