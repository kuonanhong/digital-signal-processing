 function audio_stft = mySTFT(audio_signal,window_function,step)
% Calculates Short-time Fourier transform (STFT)
% Arguments:
% audio_signal: audio signal
% window_function: window function
% step: step length in samples
           
% Number of samples and window length
number_samples = length(audio_signal);
% Calculate the total number of rows (samples)
window_length = length(window_function);

% Number of time frames. Calculate the total number of columns of audio_stft
number_times = ceil((window_length-step+number_samples)/step); 
            
% Zero-padding at the start and end to center the windows 
% in order for the configuration to be done correctly without overlapping
audio_signal = [zeros(window_length-step,1);audio_signal; ...
zeros(number_times*step-number_samples,1)];
            
% Initialize the STFT array. We define the two-dimensional array with
% window length as columns and number of time frames asrows
audio_stft = zeros(window_length,number_times);
            
% Multiplys each element of the signal with the
% hamming window and loads the results into the columns of audio_stft 
for time_index = 1:number_times
    % Window the signal
    sample_index = step*(time_index-1);
    audio_stft(:,time_index) ...
    = audio_signal(1+sample_index:window_length+sample_index).*window_function;                
end
            
% Fourier transform of audio_stft
audio_stft = fft(audio_stft);
            
end