% Sampling Frequency
Fs = 1000;
Ts = 1/Fs;
% Length of signal
L = 2001;
% Time vector for signal
t = (0:L-1)*Ts;

% White Gaussian noise (zero mean)
u = randn(L, 1);
% Implementation of dirac function at points: 0.625 & 0.650
delta = dirac(t-0.625);
delta_2 = dirac(t-0.650);
% Find Inf
idx = delta == Inf; 
% set Inf to finite value
delta(idx) = 1;
idx2 = delta_2 == Inf;
delta_2(idx2) = 1;
% Signal
x = 1.5*cos(2*pi*40*t)+1.5*cos(2*pi*100*t)+5*(delta)+5*delta_2+0.15*u.';

% Plots the signal over time
figure
plot(t, x)
title(['x(n) = 1.5 cos(2\pi40n) + 1.5cos(2\pi100n) + 0.15u(n) +'...
    ' 5(\delta(n-0.625) + \delta(n-0.650))']);
xlabel('Time(sec)');
ylabel('Amplitude');
legend('Signal');

% 4 different window lengths
winLens = [60 40 20];
nWindows = length(winLens);
k = 0;
figure
for iWinLen = winLens
    k = k+1;
    % Generate hamming window with the given length
    window = hamming(iWinLen);
    % Points that are used to calculate the stft (recommended to be power of 2) 
    nff = max(256,2^nextpow2(iWinLen));
    % Use spectrogram function to calculate the STFT of signal x with nff
    % points, sampling frequency Fs and overlap 50%. Returns a vector of 
    % cyclical frequencies F, expressed in terms of the sample rate, Fs
    [S, F, T] = spectrogram(x, window, iWinLen/2, nff, Fs);
    
    % Plots the absolute value of STFT via the function contour 
    subplot(nWindows, 1, k);
    contour(T, F, abs(S))
    ylim([0 300])
    ylabel('Frequency(Hz)')
    if (k==1)
        title('|STFT| for various Hamming window lengths')
    end
    if (k==3)
        xlabel('Time(sec)')
    end
    legend(['Window length: ',num2str(iWinLen),' Samples'])
end

% Calls wavescales (given function) to calculate vector of scales and the 
% corresponding pseudo-frequencies f
[scales, f]= wavescales('morl',Fs);
% Uses cwft (given function) to obtain the CWT coefficients using the 
% Morlet wavelet
cwtstruct = cwtft({x,1/Fs},'Scales',scales,'Wavelet','morl');
cfs = cwtstruct.cfs;
% Transpose of frequency vector f
f = f.';

% Plots the absolute value of DT-CWT in respect to frequency and time
% Uses contour function
figure
contour(t, f, abs(cfs));
title({'Wavelet Transform (Absolute value) of:'; ['x(n) = 1.5 cos(2\pi40n)'...
    ' + 1.5cos(2\pi100n) + 0.15u(n) + 5(\delta(n-0.625) + \delta(n-0.650))']});
xlabel('Time(sec)');
ylabel('Frequency(Hz)');