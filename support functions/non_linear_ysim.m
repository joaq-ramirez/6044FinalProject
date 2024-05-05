function [gamma_m] = non_linear_xsim(chi_m,w_m,w_c,n)
%non_linear_ysim Simulate EP and angular rates of spacecraft


for j = 1:2*n+1 %Iterate through chi points
    N_Rs = chi_m(:,j); % do inverse mapping from this 
    BN = EP2C(N_Rs(1:3)); 
    for i = 1:length(chi_m(:,j))
        B_Rs(i,j) = BN*N_Rs(i); 

    end
    [~, xout] = ode45(func,tspan,IC(:,j));
    chi_m(:,j) = xout(end,:)';
end
    %%%%%%
end