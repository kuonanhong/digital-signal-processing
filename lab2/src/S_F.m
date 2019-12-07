function [SF] = S_F(delta_b, P)
 % S_F calculates the region which covers the mask up until the point 
 % which undergoes the coverage
 % SF function calculates the minimum level of power that the neighboring
 % frequencies should have so that they are precieved from the human ear

 if delta_b>=-3 && delta_b<-1
    SF = 17*delta_b-0.4*P+11; 
 elseif delta_b<0
     SF = (0.4*P+6)*delta_b;
 elseif delta_b<1
     SF = -17*delta_b;
 else
     SF = (0.15*P-17)*delta_b-0.15*P;
 end
end

