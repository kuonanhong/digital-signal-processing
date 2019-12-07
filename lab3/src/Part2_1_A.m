close all, clear all

%% Part 2.1.A / Beamforming in Simulated Signals

N=7;            % number of microphones
c=340;          % speed of sound in air (m/s)
d=0.04;         % distance between the mics (m)  
theta_s=pi/4;   % angle of source

% Read audio signal input for every mic 
input(1,:)=audioread('sensor_0.wav');
input(2,:)=audioread('sensor_1.wav');
input(3,:)=audioread('sensor_2.wav');
input(4,:)=audioread('sensor_3.wav');
input(5,:)=audioread('sensor_4.wav');
input(6,:)=audioread('sensor_5.wav');
input(7,:)=audioread('sensor_6.wav');
source=audioread('source.wav');

% Sampling frequency in Hz
fs=48000;

% Angular frequency vector
w=-pi:((2*pi)/length(input(1,:))):pi-((2*pi)/length(input(1,:)));

% Fractional Delay Filter
dks = zeros(7,length(input(1,:)));
for i=1:7
    dks(i,:)=exp(j*(N-1)*w*fs*d*(cos(theta_s))/(2*c)).*exp(-j*(i-1)*w*fs*d*(cos(theta_s))/(c));
end

%Find fourier transform of every input
dftinput = zeros(7,length(input(1,:)));
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
noise = source - yfinal;
audiowrite('sim_ds.wav',yfinal,fs);

% Check if the alignment of the input signals in time is right
sound(noise,fs);

% We plot the different signals in order to compare them
k=0:(1/fs):length(source)/fs-(1/fs);
plot(k,source); % Plot the clear signal
title('Clear signal');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
k=0:(1/fs):length(yfinal)/fs-(1/fs);
plot(k,yfinal); % Plot the delay-and-sum beamformer output
title('y(t) / output of the delay-and-sum beamformer');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
k=0:(1/fs):length(input(4,:))/fs-(1/fs);
plot(k,input(4,:)); % Plot the noisy signal from the central microphone
title('Noisy signal from the central microphone');
xlabel('Time(sec)');
ylabel('Amplitude');

figure;
k=0:(1/fs):length(input(4,:))/fs-(1/fs);
plot(k,noise); % Plot the noise
title('Noise');
xlabel('Time(sec)');
ylabel('Amplitude');

% We plot the spectrogramms (PSD in logarithmic scale) 
% of the beamformer and the central microphone
% and compare them to the clear signal (source.wav)
windowlength=0.005*fs;
overlap= 0.0025 * fs;
figure;
subplot(3,1,1);
[~,F,T,P] = spectrogram(source, windowlength, overlap, nextpow2(length(source)), fs);
surf(T,F,10*log10(P),'edgecolor','none'); 
axis tight; 
view(0,90);
title('Spectogram of the clear (Source) signal');
ylabel('Hz');

subplot(3,1,2);
[~,F,T,P] = spectrogram(yfinal, windowlength, overlap, nextpow2(length(yfinal)), fs);
surf(T,F,10*log10(P),'edgecolor','none'); 
axis tight; 
view(0,90);
title('Spectogram of y(t)(output of the beamformer)');
xlabel('Time(Sec)'); 
ylabel('Hz');

subplot(3,1,3);
[S,F,T,P] = spectrogram(input(4,:), windowlength, overlap, nextpow2(length(input(4,:))), fs);
surf(T,F,10*log10(P),'edgecolor','none'); 
axis tight; 
view(0,90);
title('Spectogram of the input signal (central microphone)');
xlabel('Time (Seconds)'); 
ylabel('Hz');

% Compute the SNR of the noisy signal in the central mic and the SNR of the output
% of the beamformer
SNRcentral=snr(input(4,:), (source' - input(4,:)));
SNRbeamformer=snr(yfinal,(source - yfinal));
display(SNRcentral);
display(SNRbeamformer);
