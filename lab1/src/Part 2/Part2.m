% Audio signal averaged over its channels and sample rate in Hz
[audio_signal,sample_rate] = audioread('speech_utterance.wav');

% Window length in samples
window_duration = 0.04;
window_length = sample_rate*window_duration;

% Window function (periodic Hamming window)
window_function = hamming(window_length,'periodic');

% Step length in samples (half the window length for OLA)
step_length = window_length/2;

% Magnitude spectrogram (without the DC component and the mirrored frequencies)
audio_stft = mySTFT(audio_signal,window_function,step_length);
audio_spectrogram = abs(audio_stft(2:window_length/2+1,:));

% Recovery of original speech signal
x_istft = myISTFT(audio_stft, window_function,step_length);
sound(x_istft,sample_rate);
audiowrite('speech_utterance_rec.wav', x_istft, sample_rate);

% Length of audio signal
number_samples = length(audio_signal);
% Calculation of time vector(t) and frequency vector(Fk) for the spectogram
number_times = ceil((window_length-step_length+number_samples)/step_length);
t=(window_length:step_length:(window_length+(number_times-1)*step_length))/sample_rate;
% Fk = Wk*Fs/2pi
Wk=2*pi*(0:2:window_length-1)/window_length;
Fk=Wk*sample_rate/(2*pi);

% Spectrogram displayed in sec and Hz
figure
surf(t,Fk,audio_spectrogram)
shading interp
axis xy
colormap(jet)
title('Spectrogram of speech signal')
view(0, 90)
xlabel('Time(sec)')
ylabel('Frequency(Hz)')

% Isolation of /o/ and /a/ phonemes
% We locate the exact seconds in which the phonemes are
% and convert them in samples. In the original signal, the intervals of the
% phonemes are used to calculate their STFT and the rest of the signal is
% put to zero.

% Phoneme /o/
u=audio_signal;
for i=10501:42000
 u(i)=u(i)*0;
end
for i=43900:58720
 u(i)=u(i)*0;
end

%Phoneme /a/
y=audio_signal;
for i=1:  12300 
 y(i)=y(i)*0;
end
for i= 13590: 25300
y(i)=y(i)*0;
end
for i= 26500:58720
y(i)=y(i)*0;
end

% Perform STFT
u = mySTFT(u,window_function,step_length);
audio_spectrogram1 = abs(u(2:window_length/2+1,:));
y = mySTFT(y,window_function,step_length);
audio_spectrogram2 = abs(y(2:window_length/2+1,:));

figure
surf(t,Fk,audio_spectrogram1)
shading interp
axis xy
colormap(jet)
box on
view(0, 90)
xlabel('Time(sec)')
ylabel('Frequency(Hz)')
title('Amplitude spectrogram of the phoneme /o/')

figure
surf(t,Fk,audio_spectrogram2)
shading interp
axis xy
colormap(jet)
box on
view(0, 90)
xlabel('Time(sec)')
ylabel('Frequency(Hz)')
title('Amplitude spectrogram of the phoneme /a/')

% Time vector for the original signal
t1 = (0:number_samples-1)/sample_rate;

figure
% Plots the original signal
plot(t1,audio_signal , 'b')
grid on
xlim([0 max(t)])
xlabel('Time(sec)')
ylabel('Amplitude')
title('Original and reconstructed speech signal')
hold on

% Time vector for the reconstructed signal
number_samples_i = length(x_istft);
t2 = (0:number_samples_i-1)/sample_rate;
% plot the resynthesized signal 
plot(t2, x_istft, '-.r')
legend('Original signal', 'Reconstructed signal')
hold off
