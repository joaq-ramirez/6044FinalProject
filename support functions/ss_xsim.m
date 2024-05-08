function [chi_m] = ss_xsim(IC,c)
%non_linear_xsim Simulate EP and angular rates of spacecraft
func = @(inp)[
        -cross(inp(4:6),inp(1:3))
        c.I\(-cross(inp(4:6),c.I*inp(4:6))); 
        ];

    tstep = 1; 
    % Rsk = IC(1:3); 
    % w_sc = IC(4:6);
    
    % Rsk1 = Euler3122C(w_sc*tstep)*Rsk; 
    % Rsk1 =  -cross(w_sc,Rsk);

    k1 = tstep * func(IC(:));
    k2 = tstep * func(IC(:) + k1/2);
    k3 = tstep * func(IC(:) + k2/2);
    k4 = tstep * func(IC(:) + k3);
    chi_m= IC(:) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);

    % w_sc1 = w_sc + (c.I\(-cross(w_sc,c.I*w_sc)))*tstep; 
    % chi_m = [Rsk1; w_sc1]; 
end