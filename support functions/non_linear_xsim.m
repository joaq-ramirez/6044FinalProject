function [chi_m] = non_linear_xsim(n,IC,c)
%non_linear_xsim Simulate EP and angular rates of spacecraft

for j = 1:2*n+1 %Iterate through chi points
    func = @(t,inp)[
        inp(2); 
        -mu*inp(1)/(sqrt((inp(1)^2+inp(3)^2))^3); 
        ];
    tspan = [0 1]; %match timestep size
    [~, xout] = ode45(func,tspan,IC(:,j));
    chi_m(:,j) = xout(end,:)';
end
    %%%%%%
end