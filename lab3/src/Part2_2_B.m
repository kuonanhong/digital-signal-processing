close all, clear all

%% Postfiltering with Wiener filter

N=7;            % Number of mics
c=340;          % Speed of sound in air
d=0.04;         % Distance between the microphones
theta_s=pi/4;   % Angle of source
fs=48000;       % Sampling frequenxy

% Read input signals from every microphone
input(1,:)=audioread('sensor_0.wav');
input(2,:)=audioread('sensor_1.wav');
input(3,:)=audioread('sensor_2.wav');
input(4,:)=audioread('sensor_3.wav');
input(5,:)=audioread('sensor_4.wav');
input(6,:)=audioread('sensor_5.wav');
input(7,:)=audioread('sensor_6.wav');
source=audioread('source.wav');

% First, perform delay-sum beamforming
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
yfinal = yfinal';

% Then, we divide the signal to frames of 30ms
winlength=0.03*fs;
winoverlap=0.5*winlength;
winstep=winlength-winoverlap;
inputframe=buffer(yfinal,winlength,winoverlap); 
sourceframe=buffer(source,winlength,winoverlap);

% Short-time analysis using Hamming window
window=hamming(winlength);
numofframes=size(inputframe,2);
noisyframe=inputframe(:,3);

% We have to pass every frame from a Wiener filter
% So for every frame, we have to implement a Wiener filter
y_wiener = zeros(size(inputframe));
for j=1:size(inputframe,2)
    xt = inputframe(:,j);
    st = sourceframe(:,j );
    xt = xt .* window;
    st = st .* window;
    % Calculate PSD of every frame
    [XT,f1] = pwelch(xt,[],[],winlength,'twosided');
    [ST,f2] = pwelch(st,[],[],winlength,'twosided');
    
    Hw = (ST ./ XT); 
    xt_DFT = fft(xt);
    wiener_output = Hw.*xt_DFT;
    y_wiener(:,j) = real(ifft(wiener_output));
    
end

% Implement overlap-add method to reconstruct signal from the frames
first_half = 1:size(y_wiener,1)/2;
second_half = size(y_wiener,1)/2+1 : size(y_wiener,1);
output = y_wiener(second_half, 1:end-1) + y_wiener(first_half,2:end);
output = [y_wiener(first_half,1); output(:); y_wiener(second_half,end)];
output = output(1:length(yfinal));

audiowrite('real_mmse.wav',output,fs);
sound(output,fs);

% We plot the different signals in order to compare them
time=0:1/fs:length(yfinal)/fs-(1/fs);
figure;
plot(time,source);
title('Clear signal');
xlabel('time(sec)');
ylabel('Amplitude');

figure;
plot(time,input(4,:));
title('Noisy signal from the central microphone');
xlabel('time(sec)');
ylabel('Amplitude');

figure;
plot(time,yfinal);
title('y(t)/Output of the beamformer / Input of Wiener filter');
xlabel('time(sec)');
ylabel('Amplitude');

figure;
plot(time,real(output));
title('Output of postfilter Wiener');
xlabel('time(sec)');
ylabel('Amplitude');

% Spectrogramms of the above signals
winlength=0.005*fs;
winoverlap= 0.0025 * fs;

figure;
subplot(4,1,1);
[~,F,T,P] = spectrogram(source, winlength, winoverlap, nextpow2(length(source)), fs);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
view(0,90);
title('Spectogram of Clear signal from Source.wav');
ylabel('Hz');

subplot(4,1,2);
[~,F,T,P] = spectrogram(yfinal, winlength, winoverlap, nextpow2(length(yfinal)), fs);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
view(0,90);
title('Spectogram of y(t)/output from the beamformer / Wiener input');
ylabel('Hz');

subplot(4,1,3);
[~,F,T,P] = spectrogram(input(4,:), winlength, winoverlap, nextpow2(length(input(4,:))), fs);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
view(0,90);
title('Spectogram of input signal in central microphone');
ylabel('Hz');

subplot(4,1,4);
[~,F,T,P] = spectrogram(real(output), winlength, winoverlap, nextpow2(length(real(output))), fs);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
view(0,90);
title('Spectogram of Wiener filter output');
xlabel('Time (Seconds)'); 
ylabel('Hz');

% Compue SSNR for the input and the output of Wiener filter
winlength=0.03*fs;
winoverlap=0.5*winlength;

% Noisy frame is set to the 2nd one
noise_frame_input = yfinal(1441:2880);
frame_samples = 0.03 * 48000;
% Number of frames
frame_num = floor(length(source) / frame_samples);

sum = double(0);
num =0;
% Threshold
min = -5;

% Calculate Wiener filter input's SSNR
for i=0:(frame_num-1)
    plaisio = source( i * frame_samples +1 : (i+1) * frame_samples );
    snr_plaisiu = snr(plaisio, noise_frame_input);
    if(snr_plaisiu > 35)     % Check if more than 35
        snr_plaisiu =35;
       
    end
    if(snr_plaisiu > min)    % Check if less than threshold
    
    sum = sum + snr_plaisiu;
    num = num + 1;
      
    end
end
InputSSNR = (1/num) * sum;
display(InputSSNR);

% Calculation óu^2/noisy frame for Wiener output
noise_frame_output =  output(1441:2880);

sum1 = double(0);
num1 =0;
min1 = -5;
% Calculate Wiener filter output's SSNR
for i=0:(frame_num-1)
    plaisio = source( i * frame_samples +1 : (i+1) * frame_samples );
    snr_plaisiu = snr(plaisio, noise_frame_output);
    if(snr_plaisiu > 35)    % Check if more than 35
        snr_plaisiu =35;
       
    end
    if(snr_plaisiu > min1)  % Check if less than threshold
        
    sum1 = sum1 + snr_plaisiu;
    num1 = num1 + 1;
    end
end
OutputSSNR = (1/num1) * sum1;
display(OutputSSNR);

% Calculate the average of the SSNRs of the input signals of the
% delay-and-sum beamformer + Wienerpost-filter
InTotalSSNR=0;
for j=1:N
    sum1 = double(0);
    num1 =0;
    min3 = -5;
    input1 = input(j,:);
    noise_frame_beam =  input1(1441:2880);
    for i=0:(frame_num-1)
        plaisio = source( i * frame_samples +1 : (i+1) * frame_samples );
        snr_plaisiu = snr(plaisio, noise_frame_beam');
        if(snr_plaisiu > 35)    % Check if more than 35
           snr_plaisiu =35;
       
        end
        if(snr_plaisiu > min3)  % Check if less than threshold
        
            sum1 = sum1 + snr_plaisiu;
            num1 = num1 + 1;
        end
    end
    InTotalSSNR = InTotalSSNR + (1/num1) * sum1;
end
InTotalSSNR = InTotalSSNR/7;
display(InTotalSSNR);
