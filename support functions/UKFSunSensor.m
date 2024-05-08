function [x_ukf,P_ukf,NEES,NIS] = UKFSunSensor(x_t,y_t,tvec,c)
% load('orbitdeterm_finalproj_KFdata.mat')

xp = x_t(:,3); % initial value of x
Pp = 0.01*eye(6);  % Define covariance - 6 because only 6 DOF (quaternions are constrained to 3)
% Pp(1:3,1:3) = 0.01*eye(3);
% Rk = 0.001*eye(3); % Measurement noise - pre omega
Rk = 0.01*eye(20); % Measurement noise
Qk = 0.0001*eye(6);

% preallocate for UKF loop
n = length(Pp(:,1)); 
kappa = 0;
beta = 2; 
alpha = 0.9; % smaller values go to EKF
lambda = alpha^2 * (n+kappa)-n; 
tlast = 0;
x_ukf = zeros(6,length(tvec));
P_ukf = zeros(6,6,length(tvec));

for i = 1:length(tvec)
    %reset each loop
    % Sk = chol(Pp+Qk,'lower');
    Sk = chol(Pp,'lower');
    i
    
%% 1. dynamics prediction step from time step k->k+1
    %%%%%%%%%%%%%%%%%%%
    %part A: generate sigma points from xp
    %%%%%%%%%%%%%%%%%%%

    chik0 = xp; %Initial chi_mean
    for j = 1:n
        W = sqrt(n+lambda)*Sk(j,:)';
        chik(:,j) = xp+W; % Combining the two
    end

    for j = (n+1):2*n
        W = -sqrt(n+lambda)*Sk(j-n,:)';
        chik(:,j) = xp+W; % Combining the two
    end
    chi_comb = [chik0, chik];
    
    inp_ic = chi_comb; 
    inp_ic(1:3,:) = inp_ic(1:3,:)/ norm(inp_ic(1:3,:));
    %%%%%%%%%%%%%%%%%%%
    %part B: put sigma points through dynamic non-linear EOM
    %%%%%%%%%%%%%%%%%%%

    % [chi_m] = non_linear_xsim(inp_ic,c);
    for j = 1:2*n+1
        chi_m(:,j) = ss_xsim(inp_ic(:,j),c); 
    end

    %%%%%%%%%%%%%%%%%%%
    %part c - recombine resultants
    %%%%%%%%%%%%%%%%%%%
    for j = 1:2*n %includes 0 thus +1
        w_mi(:,j) = 1/(2*(n+lambda));
        w_ci(:,j) = w_mi(j);
    end
    w_m = [lambda/(n+lambda),w_mi];
    w_c = [lambda/(n+lambda)+1-alpha^2+beta,w_ci];

    % Averaging the Angular Velocity Portion of chi_m
    xm_p1 = sum(w_m.*chi_m,2);

    % Propogated State Estimate
    % xm_p1 = [q_bar;omega_bar]; 

    % Estimating the State Covariance    
    Wp = zeros(6,length(chi_m(1,:)));
    Pm_p1 = zeros(6);

    for j = 1:length(chi_m(1,:))
        % Differences between average state and sigma points
        Wp(:,j) = chi_m(:,j)-xm_p1;

        % State Covariance Estimate
        Pm_p1 = Pm_p1 + w_c(j)*(Wp(:,j)*Wp(:,j).');
    end
    Pm_p1 = Pm_p1 +Qk; 
%% 2. Measurement Update Step at time k+1 given observation y(k+1)
    %%%%%%%%%%%%%%%%%%%
    %part a - generate sigma pts
    %%%%%%%%%%%%%%%%%%%

    Sk_p1 = chol(Pm_p1,'lower');

    %repeat above steps for chi k+1
    chik0_p1 = xm_p1;

    for j = 1:n
        W = sqrt(n+lambda)*Sk_p1(j,:)';
        chik_p1(:,j) = xm_p1+W;
    end

    for j = (n+1):2*n
        W = -sqrt(n+lambda)*Sk_p1(j-n,:)';
        chik_p1(:,j) = xm_p1+W;
    end

    chik_p1 = [chik0_p1 chik_p1]; %combine
    chik_p1(1:3,:) = chik_p1(1:3,:)/ norm(chik_p1(1:3,:)); 
    %%%%%%%%%%%%%%%%%%%
    %part b - propagate each chi through non-linear measurment function h
    %%%%%%%%%%%%%%%%%%%
    gam_p1 = zeros(20,2*n+1); % post omega

    for j = 1:2*n+1
        [gam_p1(:,j)] = ss_ysim(chik_p1(:,j));
    end
    %%%%%%%%%%%%%%%%%%%
    %part c - get predicted measurement mean and measurement covar
    %%%%%%%%%%%%%%%%%%%   

    ym_p1(:,i) = sum(w_m.*gam_p1,2);

    Pyy_p1 = zeros(20); % post omega
    for j = 1:2*n+1
        P_iter_yy = w_c(j).*(gam_p1(:,j) - ym_p1(:,i))*(gam_p1(:,j) - ym_p1(:,i)).';
        Pyy_p1 = Pyy_p1 + P_iter_yy; 
    end
    Pyy_p1 = Pyy_p1 + Rk; 
    %%%%%%%%%%%%%%%%%%%
    %part d - get state measurement cross-covariance matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Cxy_p1 = zeros(6,20); % post omega was 6,3
    for j = 1:2*n+1
        Pxy_p1_iter = w_c(j).*(Wp(:,j)*(gam_p1(:,j)-ym_p1(:,i))');
        Cxy_p1 = Cxy_p1 + Pxy_p1_iter;
    end

    %%%%%%%%%%%%%%%%%%%
    %part e - Estimate kalman gain matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Kk_p1 = Cxy_p1*inv(Pyy_p1);

    %%%%%%%%%%%%%%%%%%%
    %part f - Perform kalman state and covariance update with observation yk+1 (nxp)
    %%%%%%%%%%%%%%%%%%%  
    
    update = Kk_p1*(y_t(i,:)' - ym_p1(:,i));

    xp_p1 = xm_p1+update;
    xp_p1(1:3) = xp_p1(1:3)/ norm(xp_p1(1:3)); 
    Pp_p1 = Pm_p1 - Kk_p1*Pyy_p1*Kk_p1'; 
    
    %Store variables
    x_ukf(:,i) = xp_p1;
    P_ukf(:,:,i) = Pp_p1;
    %Set for next iteration
    xp = xp_p1;
    Pp = Pp_p1;
    % tlast = tvec(i);
end

NEES = 1;
NIS = 1; 

figure
plot(tvec,ym_p1(:,:))
% figure
% plot(tvec,debug_update(1:3,:))
% figure
% plot(tvec,debug_update(4:6,:))

end
