function S = S_T(P)
 % S_T finds tone maskers for a given PSD array by checking for 
 % local maximum points in the various frequencies
 % If in k there is a tone masker then S(k)=1 else S(k)=0
 
 S(1) = 0;
 S(2) = 0;
 for i = 251:256
     S(i) = 0;
 end
 for i = 3:250
     [dk,num] = delta_k(i);
     for j = 1:num
      if (P(i)>P(i+1)) && (P(i)>P(i-1)) && (P(i)>(P(i+dk(j))+7)) && (P(i)>(P(i-dk(j))+7))
         S(i) = 1;
      else 
         S(i) = 0;
         break
      end
     end
 end
end

