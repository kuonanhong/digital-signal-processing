clear all;

%--------------------------------------------------------------------------
% RUNTIME < 1 minute
%--------------------------------------------------------------------------

filename = 'music_0.wav';
% reading the signal from the file and normalize 
[x, fs] = audioread (filename);
x = x/max(x);
blockLength = 512;

% Plot the original audio signal
figure
plot((1:1:length(x)),x);
title('Original audio signal (Normalised)')
xlabel('Samples')
ylabel('Amplitude')

% get size of original file to compare later
original = dir(filename);
original_size = original.bytes;

% compute how many blocks are needed and allocate the input block matrix
numBlocks = floor(length(x)/blockLength);
X = zeros(numBlocks, blockLength);

% blocking of the input audio data
for k = 1:numBlocks
    for m = 1:blockLength
        X(k,m) = x((k-1)*blockLength + m);
    end
end

% plot 60th block of samples
n=1:1:blockLength;
figure
plot(n,X(60,:));
xlim([0 512])
title('60^{th} Block of input audio signal')
xlabel('1 \leq Samples \leq 512')
ylabel('Amplitude')

% array of frequencies
FFTLength = 512;
F = [1:FFTLength/2]*(fs/FFTLength);
% initialize bark scale
b = hz2bark(F);

% initialize power spectral density arrays (P1 is temporary)
P1 = zeros(numBlocks,blockLength);
P = zeros(numBlocks,blockLength/2);
% calculate PSD array
for k = 1:numBlocks
    window = hanning(blockLength);
    window = window.';
    P1(k,:) = X(k,:).*window;
    P1(k,:) = fft(P1(k,:),blockLength);
    P1(k,:) = 10*log10((abs(P1(k,:))).^2);
    P1(k,:) = 90.302+P1(k,:);
end
% we only take half of the total frequencies because the fft's 
% result is symmetric 
for k =1:numBlocks
    P(k,:) = P1(k,1:blockLength/2);
end

% Plot power spectral density for the 60th block
figure
plot(F,P(60,:));
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('PSD for 60^{th} block')

% initialize tone and noise maskers arrays
PTM = zeros(numBlocks,blockLength/2);
PNM = zeros(numBlocks,blockLength/2);
for k = 1:numBlocks 
   % initialize and compute the boolean array S for k block
   S = zeros(blockLength/2); 
   P2 =P(k,:); 
   S = S_T(P2);
   % compute the power for every tone masker
   for i = 1:256
       if S(i) == 0
           PTM(k,i) = 0;
       else
           PTM(k,i) = 10*log10(10^(0.1*P(k,i))+10^(0.1*P(k,i-1))+10^(0.1*P(k,i+1)));
       end
   end
   % compute noise maskers for k block
   PNM(k,:) = findNoiseMaskers(P(k,:), PTM(k,:), b);
end

% calculate the Absolute Threshold of Hearing
Tq = hz2dBSPL(F);

PTM_temp = PTM(60,:);
PNM_temp = PNM(60,:);
PTM_temp(PTM_temp == 0) = NaN;
PNM_temp(PNM_temp == 0) = NaN;
% Plot tone and noise maskers for the 60th block in barks
figure
ax1 = subplot(2,1,1);
plot(ax1,b,P(60,:));
hold on
plot(ax1,b,PTM_temp,'*');
hold on
plot(ax1,b,Tq,'--');
hold off
ylabel(ax1,'Magnitude')
ylim(ax1,[-50 250])
title(ax1,'Tone Maskers for 60^{th} block')
legend(ax1,'Original PSD','Tone Maskers','ATH')

ax2 = subplot(2,1,2);
plot(ax2,b,P(60,:));
hold on
plot(ax2,b,PNM_temp,'o');
hold on
plot(ax2,b,Tq,'--');
hold off
ylabel(ax2,'Magnitude')
ylim(ax2,[-50 250])
title(ax2,'Noise Maskers for 60^{th} block')
legend(ax2,'Original PSD','Noise Maskers','ATH')
xlabel('Frequency (Bark)')

% initialize Global Masking Threshold array
Tg = zeros(numBlocks,blockLength/2);
for k = 1:numBlocks
    % reorganization of tone and noise maskers via given function
    [PTM(k,:), PNM(k,:)] = checkMaskers(PTM(k,:), PNM(k,:) , Tq, b);
    
    % initialize Individual Masking Threshold arrays
    TTM = zeros(blockLength/2,blockLength/2);
    TNM = zeros(blockLength/2,blockLength/2);
    for j = 1:256
       % if we find tone masker
       if PTM(k,j)>0  
           % then for every discrete frequency 
            for i = 1:256              
                % if we are in the 12-Bark neighborhood
                if b(i)>=b(j)-3 && b(i)<=b(j)+8
                    % find the distance between the frequencies in Bark
                    delta_b = b(i)-b(j); 
                    SF = S_F(delta_b,PTM(k,j));
                   % calculate Individual Masking Threshold for tone masker
                    TTM(i,j) = PTM(k,j) - 0.275*b(j) + SF - 6.025;
                    % add it to the Global Masking Threshold of frequency i
                    Tg(k,i) = Tg(k,i) + 10^(0.1*TTM(i,j));
                else
                    TTM(i,j) = 0;
                    Tg(k,i) = Tg(k,i) + 10^(0.1*TTM(i,j));
                end
            end
       end 
       % if we find noise masker
       if PNM(k,j)>0
           % then for every discrete frequency 
            for i = 1:256
                % if we are in the 12-Bark neighborhood
                if b(i)>=b(j)-3 && b(i)<=b(j)+8
                    % find the distance between the frequencies in Bark
                    delta_b = b(i)-b(j);
                    SF = S_F(delta_b,PNM(k,j));
                    % calculate Individual Masking Threshold for noise masker
                    TNM(i,j) = PNM(k,j) - 0.175*b(j) + SF - 2.025;
                    % add it to the Global Masking Threshold of frequency i
                    Tg(k,i) = Tg(k,i) + 10^(0.1*TNM(i,j)); 
                else
                    TNM(i,j) = 0;
                    Tg(k,i) = Tg(k,i) + 10^(0.1*TNM(i,j)); 
                end
            end
       end
        Tg(k,j) = Tg(k,j) + 10^(0.1*Tq(j));
    end
        % calculate the total Global Masking Threshold of block k
        Tg(k,:) = 10*log10(Tg(k,:));
end

check_PTM_temp = PTM(60,:);
check_PNM_temp = PNM(60,:);
check_PTM_temp(check_PTM_temp == 0) = NaN;
check_PNM_temp(check_PNM_temp == 0) = NaN;
% Plot checked tone and noise maskers for the 60th block in barks
figure
ax1 = subplot(2,1,1);
plot(ax1,b,P(60,:));
hold on
plot(ax1,b,check_PTM_temp,'*');
hold on
plot(ax1,b,Tq,'--');
hold off
ylabel(ax1,'Magnitude')
ylim(ax1,[-50 250])
title(ax1,'Checked Tone Maskers for 60^{th} block')
legend(ax1,'Original PSD','Checked Tone Maskers','ATH')

ax2 = subplot(2,1,2);
plot(ax2,b,P(60,:));
hold on
plot(ax2,b,check_PNM_temp,'o');
hold on
plot(ax2,b,Tq,'--');
hold off
ylabel(ax2,'Magnitude')
ylim(ax2,[-50 250])
title(ax2,'Checked Noise Maskers for 60^{th} block')
legend(ax2,'Original PSD','Checked Noise Maskers','ATH')
xlabel('Frequency (Bark)')

% plot Global Masking Threshold for the 60th block
figure
plot(b,P(60,:));
hold on
plot(b,check_PTM_temp,'*');
hold on
plot(b,check_PNM_temp,'o');
hold on
plot(b,Tq,'--');
hold on
plot(b,Tg(60,:));
hold off
ylabel('Magnitude')
ylim([-50 250])
title('Global Masking Threshold for 60^{th} block')
legend('Original PSD','Tones','Noise','ATH','Overall Threshold')
xlabel('Frequency (Bark)')

% number of filters
M = 32;
% length of every filter
L = 2*M;
% initialize filter arrays
F_analysis = zeros(M,L);
F_synthesis = zeros(M,L);
% calculate analysis filters
for k = 0:M-1
    for n = 0:L-1
        F_analysis(k+1,n+1) = sin((n+1/2)*pi/(2*M))*sqrt(2/M)*cos((2*n+M+1)*(2*k+1)*pi/(4*M));
    end
end
% calculate synthesis filters using the analysis filters
for k = 1:M
   for i = 0:L-1
       F_synthesis(k,i+1) = F_analysis(k,2*M-i);
   end
end

% initialize sum of all bits used
bk_sum = 0;
% Call Part2 function to get the final reconstructed signal in block form
% ------------------------------------------------------------------------
% PLEASE SELECT THE DESIRED QUANTIZATION METHOD
% FOR ADAPTIVE QUANTIZER SELECT 1 FOR NON-ADAPTIVE QUANTIZER SELECT 0
 adap_or_not = 1;
% ------------------------------------------------------------------------
for k = 1:numBlocks
    [X_Rebuild(k,:),bk_sum_temp] = Part2(X(k,:),F_analysis,F_synthesis,Tg(k,:),adap_or_not);
    bk_sum = bk_sum + bk_sum_temp;
end

% the final reconstruction of the output signal is performed using the
% overlap_add function
s = overlap_add(X_Rebuild);
% use only the useful part of the vector
s = nonzeros(s);

% create compressed music file
% For Adaptive Quantizer the filename is music_1
% For Non-Adaptive Quantizer the filename is music_2
if adap_or_not == 1
 filename_new = 'music_1.wav';
else
 filename_new = 'music_2.wav';
 % We normilize this value because it's below -1 and audiowrite produces
 % a Warning. If we do not, MATLAB will remove the sample from the signal
 s(191491) = -0.999999999;
end
 audiowrite(filename_new,s,fs);
% get size of compressed signal
fin = dir(filename_new);
fin_size = fin.bytes;

% calculate the percentage of compression
compPercent = bk_sum/M/numBlocks/16;

% Calculate the Mean Square Error of the original and the compressed signal
% using the immse built-in function.
% We advance the signal by a factor of 2*M=64 samples due to the delay
% the synthesis filters add
% We use 1:650273 samples of the original so that the MSE can be calculated
% correctly (it needs both signals to be of equal length)
mse = immse(s(64:650336),x(1:650273));

% Plot the reconstructed-compressed signal
figure
plot((1:1:length(s)),s);
title('Compressed audio signal (Normalised)')
xlabel('Samples')
ylabel('Amplitude')

% Plot the error between the original and the compressed signal
figure
plot((1:1:650273),-s(64:650336)+x(1:650273));
title('Original-Compressed Signal Error')
xlabel('Samples')
ylabel('Amplitude')
if adap_or_not == 0
 ylim([-0.03 0.03])
else
 ylim([-0.008 0.008])
end

% -------------------------------------------------------------------------
% Uncomment to listen the compressed signal
% sound(s,fs);
% -------------------------------------------------------------------------