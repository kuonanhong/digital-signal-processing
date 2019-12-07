clear all, close all

%% Part 2.2.A / Beamforming in Real Signals

N=7;                    % Number of mics
c=340;                  % Speed of sound in air
d=0.04;                 % Distance between microphones
theta_s=pi/4;           % Angle of source
fs=48000;               % Sampling frequency

% Read input signals from every microphone
input(1,:)=audioread('sensor_0.wav');
input(2,:)=audioread('sensor_1.wav');
input(3,:)=audioread('sensor_2.wav');
input(4,:)=audioread('sensor_3.wav');
input(5,:)=audioread('sensor_4.wav');
input(6,:)=audioread('sensor_5.wav');
input(7,:)=audioread('sensor_6.wav');
source=audioread('source.wav');

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
audiowrite('real_ds.wav',yfinal,fs);

% Check if the alignment of the input signals in time is right
%sound(noise,fs);

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

% Central microphone
central_mic = input(4,:);
% Noisy frame is set to the 2nd one
noise_frame_central = central_mic(1441:2880);
frame_samples = 0.03 * 48000;
% Number of frames
frame_num = floor(length(source) / frame_samples);

sum = double(0);
num =0;
% Threshold
min = 0;

% Calculate central microphone's SSNR
for i=0:(frame_num-1)
    plaisio = source( i * frame_samples +1 : (i+1) * frame_samples );
    snr_plaisiu = snr(plaisio, noise_frame_central');
    if(snr_plaisiu > 35)     % Check if more than 35
        snr_plaisiu =35;
       
    end
    if(snr_plaisiu > min)    % Check if less than threshold
    
    sum = sum + snr_plaisiu;
    num = num + 1;
      
    end
end
SSNR_central = (1/num) * sum;
display(SSNR_central);

% Calculation óu^2/noisy frame for beamformer output
noise_frame_beam =  yfinal(1441:2880);

sum1 = double(0);
num1 =0;
min = 0;
% Calculate beamformer's output SSNR
for i=0:(frame_num-1)
    plaisio = source( i * frame_samples +1 : (i+1) * frame_samples );
    snr_plaisiu = snr(plaisio, noise_frame_beam);
    if(snr_plaisiu > 35)  % Check if more than 35
        snr_plaisiu =35;
       
    end
    if(snr_plaisiu > min) % Check if less than threshold
        
    sum1 = sum1 + snr_plaisiu;
    num1 = num1 + 1;
    end
end
SSNR_beamformer = (1/num1) * sum1;

display(SSNR_beamformer);
