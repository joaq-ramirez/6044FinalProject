function [B3] = EP_Add(B1,B2)
    
% EP_Add	
%
%	B3 is the resulting attitude after performing rotation B1
%	followed by rotation B2. B3, B2, and B1 are all Euler 
%   Parameter attitude descriptions.
%
    
    B3 = [B2(1) -B2(2) -B2(3) -B2(4);
          B2(2) B2(1) B2(4) -B2(3);
          B2(3) -B2(4) B2(1) B2(2);
          B2(4) B2(3) -B2(2) B2(1)]*B1;

    B3 = B3/norm(B3);
end