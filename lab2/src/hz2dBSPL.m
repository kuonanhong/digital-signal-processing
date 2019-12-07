function [Tq] = hz2dBSPL(hz)
 % HZ2DBSPL calculates the Absolute Threshold of Hearing
 % changes frequency in Hz to frequency in dB SPL
 % via the following equation
Tq = 3.64*(hz/1000).^(-0.8)-6.5*exp(-0.6*(hz/1000-3.3).^2)+10^(-3)*(hz/1000).^4;
end

