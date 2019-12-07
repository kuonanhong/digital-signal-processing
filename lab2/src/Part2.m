function [X_Rebuild,bk_sum_temp] = Part2(X,F_analysis,F_synthesis,Tg,adap_or_not)
 % Part2 filters the given block of the original signal using the
 % filterbank specified by F_analysis and F_synthesis filters
 % The output is the reconstructed signal and the number of bits used
 % You can SELECT ONLY ONE method out of two: 
 % 1)Adaptive Unifrom Quantizer
 % 2)Non-Adaptive Uniform Quantizer
 
 % number of filters
 M = 32;
 % number of levels of original signal
 R = 2^16;
 % output of synthesis filter
 X_new = zeros(M,639);
 % output of filterbank (sum of all X_new)
 X_Rebuild = zeros(1,639);
 % sum of all bits used
 bk_sum_temp = 0;
 
 for k = 1:M
    % convolution of original signal and analysis filter
    v = conv(X,F_analysis(k,:));
    % downsampling by a factor of M
    y = downsample(v,M);
    
    % split Tg(1:256) into regions of 8 samples (32*8=256)
    % so that every filter uses exaclty one region of Tg to calculate bk
    indexes = (8*(k-1)+1):8*k; 
    
    % Perform quantization using two methods(SELECT ONLY ONE)
    % Quantization parameters
    if adap_or_not == 1
     % SELECT FOR ADAPTIVE UNIFORM QUANTIZER 
     bk = ceil(log2(R./min(Tg(indexes))) - 1);
     range = max(y) -min(y);
    else
     % SELECT FOR NON-ADAPTIVE UNIFORM QUANTIZER
     bk = 8;
     range = 2;
    end
    bk_sum_temp = bk_sum_temp + bk;
    
    % quantization step
    delta = range/2^bk;
    % initialize quantized signal
    y_new = zeros(1,length(y));
    % Perform quantization using the parameters specified above
    q = range/(2^bk-1);
    h = 2*floor(y./q)+1;
    y_new = h.*delta/2;
    
    % upsampling by a factor of M
    wo = upsample(y_new,M);
    % convolution of wo and synthesis filter
    X_new(k,:) = conv(F_synthesis(k,:),wo);
    % calculate the output of filterbank
    X_Rebuild = X_Rebuild + X_new(k,:); 
 end
end
