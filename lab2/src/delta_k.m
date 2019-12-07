function [dk,num] = delta_k( k )
 % DELTA_K returns a set of integers for every k given
 % , where k is a discrete frequency
 % the return value is either 2, [2,3] or [2,6]
 if (k>2 && k<63)
     dk = 2;
     num  = 1;
 elseif (k>62 && k<127)
     dk(1) = 2;
     dk(2) = 3;
     num = 2;
 else
     dk(1) = 2;
     dk(2) = 3;
     dk(3) = 4;
     dk(4) = 5;
     dk(5) = 6;
     num = 5;
 end
end

