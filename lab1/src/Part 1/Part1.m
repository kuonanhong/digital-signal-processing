% Computation of ST-ZCR and STE of speech/music signals.
% Functions required: zerocross, sgn, winconv.

% Read the speech signal and save it's sampling frequency in Fs
% [x,Fs] = audioread('speech_utterance.wav');
 [x,Fs] = audioread('music.wav');

% Transpose of signal
x = x.';
% Signal length
N = length(x);
n = 0:N-1;
% Time vector
ts = n*(1/Fs); 

% Define the window
wintype = 'hamming';
win_duration = 0.02;
winlen = win_duration*Fs+1;
winamp = [0.5,1]*(1/winlen);

% Calculate the zero-crossing rate
zc = zerocross(x,wintype,winamp(1),winlen);

% Calculate the short-time energy
E = energy(x,wintype,winamp(2),winlen);

% Time index for the ST-ZCR and STE after delay compensation
% STE/ZTE are delayed due to lowpass filtering. This delay is
% compensated for the graph.
out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
t = (out-(winlen-1)/2)*(1/Fs);

% Voiced-Unvoiced Speech for speech_utterance.wav
% Corresponds to 'OLA'
% voiced_1 = x(Fs*0.6:Fs*0.82);
% sound(voiced_1,Fs);
% Corresponds to 'A H'
% voiced_2 = x(Fs:Fs*1.15);
% sound(voiced_2,Fs);
% Corresponds to 'AN H'
% voiced_3 = x(Fs*1.24:Fs*1.32);
% sound(voiced_3,Fs);
% Corresponds to 'ALO'
% voiced_4 = x(Fs*2.5:Fs*2.8);
% sound(voiced_4,Fs);

% Corresponds to 'FT'
% unvoiced_1 = x(Fs*0.85:Fs);
% sound(unvoiced_1,Fs);
% Corresponds to 'T'
% unvoiced_2 = x(Fs*1.23:Fs*1.34);
% sound(unvoiced_2,Fs);
% Corresponds to 'M'
% unvoiced_3 = x(Fs*1.55:Fs*1.67);
% sound(unvoiced_3,Fs);
% Corresponds to 'ST'
% unvoiced_4 = x(Fs*2.1:Fs*2.25);
% sound(unvoiced_4,Fs);

figure;
plot(ts,x);
title('Speech/Music Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
legend('Original Signal');

figure
plot(t,zc(out),'r'); 
xlabel('Time(sec)');
title(['Short-time Zero Crossing Rate/Window length=' num2str(win_duration) 'msec']);
legend('STZCR');

figure;
plot(t,E(out),'r'); 
xlabel('Time(sec)');
title(['Short-time Energy/Window length=' num2str(win_duration) 'msec']);
legend('STE');