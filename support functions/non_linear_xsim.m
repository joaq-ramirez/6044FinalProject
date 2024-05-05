function [chi_m] = non_linear_xsim(IC,c)
%non_linear_xsim Simulate EP and angular rates of spacecraft
c.tstep = 1; 
func = @(inp)[
        0.5*BmatQuat(inp(1:4))*[inp(5:7)]; 
        c.I\cross(-inp(5:7),c.I*inp(5:7)); 
        ];
tspan = [0 c.tstep]; %match timestep size

n = length(IC(1,:));
chi_m = zeros(7,n);

for j = 1:n % Iterate through chi points
%     [~, xout] = ode45(func,tspan,IC(:,j));

    k1 = c.tstep * func(IC(:,j));
    k2 = c.tstep * func(IC(:,j) + k1/2);
    k3 = c.tstep * func(IC(:,j) + k2/2);
    k4 = c.tstep * func(IC(:,j) + k3);
    xout= IC(:,j) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);

    % Re-Normalize the Quaternion
    xout(1:4) = xout(1:4)/norm(xout(1:4));

    chi_m(:,j) = xout.'; % prediction
end


    %%%%%%
end