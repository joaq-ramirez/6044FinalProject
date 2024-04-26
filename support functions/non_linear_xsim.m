function [chi_m] = non_linear_xsim(n,IC,c)
%non_linear_xsim Simulate EP and angular rates of spacecraft
tstep = 1; 
func = @(inp)[
        0.5*BmatQuat(inp(1:4))*[inp(5:7)]; 
        c.I\cross(-inp(5:7),c.I*inp(5:7)); 
        ];
tspan = [0 1]; %match timestep size

for j = 1:2*n+1 %Iterate through chi points
    % [~, xout] = ode45(func,tspan,IC(:,j));

    k1 = tstep * func(IC(:,j));
    k2 = tstep * func(IC(:,j) + k1/2);
    k3 = tstep * func(IC(:,j) + k2/2);
    k4 = tstep * func(IC(:,j) + k3);
    xout(:,j+1) = IC(:,j) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);

    chi_m(:,j) = xout(:,j)'; % prediction
end
    %%%%%%
end