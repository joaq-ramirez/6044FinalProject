function [chi_m] = non_linear_xsim(n,IC,c)
%non_linear_xsim Simulate EP and angular rates of spacecraft

func = @(t,inp)[
        0.5*BmatQuat(inp(1:4))*[0 ; inp(5:7)]; 
        c.I\(-inp(5:7)'*c.I*inp(5:7)); 
        ];
tspan = [0 1]; %match timestep size

for j = 1:2*n+1 %Iterate through chi points
    [~, xout] = ode45(func,tspan,IC(:,j));
    chi_m(:,j) = xout(end,:)';
end
    %%%%%%
end