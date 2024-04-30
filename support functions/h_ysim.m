function [gamma_m] = h_ysim(chik_p1)
%non_linear_ysim Map x state to y measuremnts of the sun vector

N_Rs = [1; 0; 0]; %% double check that this is right

beta_BN = chik_p1(1:4);
BN = EP2C(beta_BN);

gamma_m = BN*N_Rs; 

end