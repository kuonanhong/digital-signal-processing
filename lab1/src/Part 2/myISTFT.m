 function audio_signal = myISTFT(stft,window_function,step)
% Calculates Inverse Short-Time Fourier transform (STFT)   
% Arguments:
% audio_stft: STFT of audio signal
% window_function: window function 
% step: step length in samples
% audio_signal: reconstructed audio signal 

% Window length in samples and number of time frames
[window_length,number_times] = size(stft);

% Number of samples for the reconstructed signal using the total number of 
% time frames, step length and window length (inverse process of number_times)
number_samples = (number_times-1)*step+window_length;
            
% Initialize the reconstructed signal. We define the one-dimensional array
% that will store the reconstructed signal
 audio_signal = zeros(number_samples,1);

% Inverse Fourier transform of the stft. We use real function to
% make sure we end up with real values
i_stft = real(ifft(stft));
 
% Implementation of OLA technique with 50% overlay. we add the inverse transform
for time_index = 1:number_times
    % Constant overlap-addition
    sample_index = step*(time_index-1);
    % Constant addition of the inverse trasform
    audio_signal(1+sample_index:window_length+sample_index) ...
    = audio_signal(1+sample_index:window_length+sample_index)+i_stft(:,time_index);
end

% Remove the zero-padding at the start and at the end which we had added 
% to the beginning of the configuration
audio_signal = audio_signal(window_length-step+1:number_samples-(window_length-step));

% Division with the sum of all the windows
audio_signal = audio_signal/sum(window_function(1:step:window_length));

 end