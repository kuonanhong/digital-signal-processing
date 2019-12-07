function barks = hz2bark(hz)
 % hz2bark changes frequency in Hz to frequency in Barks
 % via the following equation
 barks=13*atan(.00076*hz)+3.5*atan((hz/7500).^2); 
end

