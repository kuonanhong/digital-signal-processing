% Sampling Frequency
Fs = 1000;
Ts = 1/Fs;
% Length of signal
L = 2001;
% Time vector for signal
t = (0:L-1)*Ts;

% White Gaussian noise (zero mean)
u = randn(L, 1);
% Signal
x = 1.5*cos(2*pi*80*t)+2.5*sin(2*pi*150*t)+0.15*u.';

% Plots the signal over time
figure
plot(t, x);
title('x(n) = 1.5 cos(2\pi80n) + 2.5sin(2\pi150n) + 0.15u(n)');
xlabel('Time(sec)');
ylabel('Amplitude');
legend('Signal');

% Window length
winLen = 40;
% Generate hamming window
window = hamming(winLen);
% Points that are used to calculate the stft (recommended to be power of 2) 
nff = max(256,2^nextpow2(winLen));
% Use spectrogram function to calculate the STFT of signal x with nff
% points, sampling frequency Fs and overlap 50%. Returns a vector of 
% cyclical frequencies F, expressed in terms of the sample rate, Fs
[S, F, T] = spectrogram(x, window, winLen/2, nff, Fs);
plot_data_1 = abs(S);

% Plots the absolute value of STFT via the function surf (can be viewd as a
% 2-D or 3-D graph
figure
h_1 = surf(T, F, plot_data_1);
set(h_1, 'edgecolor', 'none')
az = 0;
el = 90;
view(az, el);
colorbar;
xlabel('Time(sec)');
ylabel('Frequency(Hz)');
zlabel('Magnitude');
title({'Spectrogram of:'; ' x(n) = 1.5 cos(2\pi80n) + 2.5sin(2\pi150n) + 0.15u(n)'});

% Calls wavescales (given function) to calculate vector of scales and the 
% corresponding pseudo-frequencies f
[scales, f]= wavescales('morl',Fs);
% Uses cwft (given function) to obtain the CWT coefficients using the 
% Morlet wavelet
cwtstruct = cwtft({x,1/Fs},'Scales',scales,'Wavelet','morl');
cfs = cwtstruct.cfs;
% Transpose of frequency vector f
f = f.';
plot_data_2 = abs(cfs);

% Plots the absolute value of DT-CWT in respect to frequency and time
% Uses surf function (can be viewed as 2-D or 3-D plot)
figure
h_2 = surf(t, f, plot_data_2);
set(h_2, 'edgecolor', 'none')
colorbar;
az = 0;
el = 90;
view(az, el);
xlabel('Time(sec)');
ylabel('Frequency(Hz)');
zlabel('Magnitude');
title({'Wavelet Transform of:'; ' x(n) = 1.5 cos(2\pi80n) + 2.5sin(2\pi150n) + 0.15u(n)'});