function [s] = overlap_add(X)
 % overlap_add calculates output signal s using the OverLap-Add method 
 % Because of the convolution used in Part2 the signal X has more samples
 % than the original, so we make sure every block starts at the 513 sample
 % of the previous one (the original signal has blocks of 512 samples)
 % The in-between samples are of-course added together
 
 % initialization of output signal
 s = zeros(700000,1);
 X = X.';

 % the output signal is constructed using the blocks of the array X
 % using the method explained in the beginning
 j =1;
 for i = 1:1270
     len = j+608-1;
     s(j:len) = s(j:len) + X(1:608,i);
     j = len-95;
 end
end

