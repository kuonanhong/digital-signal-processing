clear all, close all

%% Part 2.1.B / Single Channel Speech Enhancement Using Wiener Filtering

N=7;                            % number of microphones
c=340;                          % speed of sound in air
d=0.04;
theta_s=pi/4;                   % angle of source
fs=48000;                       % sampling frequency
fstart = 0.36;                  % starting/ending point of frame
fend = 0.39;
frame_length = fend - fstart;  % 30msec window
frame_samples = fs*frame_length;
frame_num = fstart / (fend-fstart);
winlength=0.005*fs;

% Read input signals from the microphones
input(1,:)=audioread('sensor_0.wav');
input(2,:)=audioread('sensor_1.wav');
input(3,:)=audioread('sensor_2.wav');
input(4,:)=audioread('sensor_3.wav');
input(5,:)=audioread('sensor_4.wav');
input(6,:)=audioread('sensor_5.wav');
input(7,:)=audioread('sensor_6.wav');
source=audioread('source.wav');

% Extract the frame samples from the central microphone's and the source's
% signal
central_mic = input(4,:)';
xt = central_mic( fstart*fs : fend*fs );
st = source( fstart*fs : fend*fs ); 

% Calculate the PSD of noise u(t)
ut = xt - st;
[UT] = pwelch(ut,[],[],1441,'twosided');

% IIR Wiener response using PSD of x(t) and s(t) signals
[XT] = pwelch(xt,[],[],1441,'twosided');
[ST] = pwelch(st,[],[],1441,'twosided');
Hw = ST./XT;   

% Plot the filter response in logarithmic scale for frequencies of [0,8]kHz
f = linspace(0,8000,240);
figure
plot(f, (10*log10(Hw(1:240))));
title('IIR Wiener Filter Response');
xlabel('Frequency(Hz)'); 
ylabel('dB');

% Calculate the speech distortion index
nw = (abs(1-Hw)).^2;

% And plot it in respect to frequency
figure
plot(f, 10*log10(nw(1:240))) ;
title('Speech distortion index');
xlabel('Frequency(Hz)'); 
ylabel('dB');

% Calculate the Fourier transform of the noisy input signal x(t)
Xt = fft(xt);

% Calculate the output of the Wiener filter
wiener_output = Hw .* Xt;
xt_filtered = ifft(wiener_output);

% Calculate the PSD of the filtered signal
[XT_filtered] = XT.*((abs(Hw)).^2);

% Plot the PSD's of the input signal (with and without noise) and the output
% of the Wiener filter
figure
plot(f, 10*log10(ST(1:240)), 'k', f, 10*log10(XT(1:240)), 'r', f, 10*log10(XT_filtered(1:240)), 'b', f, 10*log10(UT(1:240)), 'g');
legend('s(t)', 'x(t)', 'y(t)', 'u(t)');
title('Input and Output signals / Wiener Filtering');
xlabel('Frequency(Hz)'); 
ylabel('dB');

% Compare SNR of input and output signals
yt2 = real(xt_filtered);
snr_xt = snr(xt, ut);
SNR_xt = snr_xt;
ut2 = yt2-st;
snr_yt = snr(yt2, ut2);
SNR_Wiener = snr_yt;

display(SNR_xt);
display(SNR_Wiener);

% Émplement multi-channel method to compare with the signle-channel method

% Angular frequency vector
w=-pi:((2*pi)/length(input(4,:))):pi-((2*pi)/length(input(4,:)));

% Fractional Delay Filter
dks = zeros(7,length(input(4,:)));
for i=1:7
    dks(i,:)=exp(j*(N-1)*w*fs*d*(cos(theta_s))/(2*c)).*exp(-j*(i-1)*w*fs*d*(cos(theta_s))/(c));
end

%Find fourier transform of every input
dftinput = zeros(7,length(input(4,:)));
for i=1:7
    dftinput(i,:)=fft(input(i,:));
end

% Pass every input from the delay filters
Y = zeros(7,length(input(1,:)));
y = zeros(7,length(input(1,:)));
yreal = zeros(7,length(input(1,:)));
for i=1:7
    Y(i,:)=(dks(i,:)) .*(dftinput(i,:));        
    y(i,:)=ifft(Y(i,:));
    yreal(i,:)=(y(i,:)+conj(y(i,:)))/2; %take only the real part of the signal
end

% Sum all the filtered inputs to get the final output
yfinal=(yreal(1,:)+yreal(2,:)+yreal(3,:)+yreal(4,:)+yreal(5,:)+yreal(6,:)+yreal(7,:))/N;

% Take the difference between the beam formed signal and the clear signal 
yfinal = yfinal';

% Get the SNR of the frame being processed
xt2 = yfinal(fs*fstart:fs*fend);
ut3 = xt2 - st;
[UT3] = pwelch(ut3,[],[],1441,'twosided');
[XT2] = pwelch(xt2,[],[],1441,'twosided');
snr_xt2 = snr(xt2, ut3);
SNR_Beamformer = snr_xt2;
display(SNR_Beamformer);

% And display the difference between the two methods
Difference = snr_xt2 - snr_yt;
display(Difference);

% Plot the PSD's of the input frame (with and without noise) and the output
% frame of the Wiener filter
figure
plot(f, 10*log10(ST(1:240)), 'k', f, 10*log10(XT(1:240)), 'r', f, 10*log10(XT2(1:240)), 'b', f, 10*log10(XT_filtered(1:240)), 'g');
legend('s(t)', 'x(t)', 'Beamformer output', 'Wiener output');
title('Input signals / Wiener output / Beamformer output');
xlabel('Frequency(Hz)'); 
ylabel('dB');
