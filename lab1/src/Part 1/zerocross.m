function zc = zerocross(x,wintype,winamp,winlen)
% Zero Crossing Rate computation.
% Computes zero crossing rate of the sequence X. 
% wintype defines the window type.
% winamp sets the amplitude of the window 
% wlen is the length of the window

% Generate x[n] and x[n-1]
x1 = x;
x2 = [0, x(1:end-1)];

% Generate the first difference
firstDiff = sgn(x1)-sgn(x2);

% Calculate the magnitude
absFirstDiff = abs(firstDiff);

% Lowpass filtering with window
zc = winconv(absFirstDiff,wintype,winamp,winlen);