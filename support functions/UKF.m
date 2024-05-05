function [x_ukf,P_ukf,NEES,NIS] = UKF(x_t,y_t,tvec,c)
% load('orbitdeterm_finalproj_KFdata.mat')

xp = x_t(:,1); % initial value of x
Pp = 0.01*eye(6);  % Define covariance - 6 because only 6 DOF (quaternions are constrained to 3)
% Pp(1:3,1:3) = 0.01*eye(3);
Rk = 0.001*eye(3); % Measurement noise
Qk = 1e-3*eye(6); % Will need to adjust later

% preallocate for UKF loop
n = length(Pp(:,1)); 
% kappa = 0;
% beta = 0; 
% alpha = sqrt(2); % smaller values go to EKF
% lambda = alpha^2 * (n+kappa)-n; 
% tlast = 0;
x_ukf = zeros(7,length(tvec));
P_ukf = zeros(6,6,length(tvec));

for i = 1:length(tvec)
    %reset each loop
    Sk = chol(Pp+Qk,'lower');
    i;
    
%% 1. dynamics prediction step from time step k->k+1
    %%%%%%%%%%%%%%%%%%%
    %part A: generate sigma points from xp
    %%%%%%%%%%%%%%%%%%%

    chik0 = xp; %Initial chi_mean
    for j = 1:n
        % W = sqrt(n+lambda)*Sk(j,:)';
        W = sqrt(2*n)*Sk(j,:)';

        theta = norm(W(1:3)); % Principal Rotation of perturbations on attitude
        if theta < 1e-8
            q = [1 0 0 0].';
        else
            e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
            q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
            q = q/norm(q);
        end
        chik(:,j) = [EP_Add(xp(1:4),q); xp(5:7)+W(4:6)]; % Combining the two
    end

    for j = (n+1):2*n
        % W = -sqrt(n+lambda)*Sk(j-n,:)';
        W = -sqrt(2*n)*Sk(j-n,:)';
        theta = norm(W(1:3)); % Principal Rotation of perturbations on attitude
        if theta < 1e-8
            q = [1 0 0 0].';
        else
            e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
            q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
            q = q/norm(q);
        end
        chik(:,j) = [EP_Add(xp(1:4),q); xp(5:7)+W(4:6)]; % Combining the two
    end
    chi_comb = [chik0, chik];
    
    inp_ic = chi_comb; 
    %%%%%%%%%%%%%%%%%%%
    %part B: put sigma points through dynamic non-linear EOM
    %%%%%%%%%%%%%%%%%%%

    [chi_m] = non_linear_xsim(inp_ic,c);

    %%%%%%%%%%%%%%%%%%%
    %part c - recombine resultants
    %%%%%%%%%%%%%%%%%%%
%     for j = 1:2*n %includes 0 thus +1
%         w_mi(:,j) = 1/(2*(n+lambda));
%         w_ci(:,j) = w_mi(j);
%     end
%     w_m = [lambda/(n+lambda),w_mi];
%     w_c = [lambda/(n+lambda)+1-alpha^2+beta,w_ci];
%     
%     % add quaternions for xm_p1 sum
%     xmprev = eye(3); 
%     for j = 1:2*n+1
%         xmcurr = EP2C(chi_m(1:4,j)); % what is the perturbed quaternion???
%         xm_p1_mat = xmcurr*xmprev;  % separate chi_m into quaternions and angular rates then sum with weights ---
%         xmprev = EP2C(chi_m(1:4,j));
%     end
%     xm_p1 = sum(w_m.*chi_m,2); % sum along dim 2
%     xm_p1(1:4) = C2EP(xm_p1_mat); 

    % Averaging the Quaternion portion of chi_m and recording error vectors
    % (for the Covariance Calcualtions)
    [q_bar, E] = Q_mean(chi_m(1:4,:));

    % Averaging the Angular Velocity Portion of chi_m
    w_bar = 1/(2*n+1)*sum(chi_m(5:7,:),2);

    % Propogated State Estimate
    xm_p1 = [q_bar;w_bar]; 

    % Estimating the State Covariance    
    Wp = zeros(6,length(chi_m(1,:)));
    Pm_p1 = zeros(6);

    for j = 1:length(chi_m(1,:))
        % Differences between average state and sigma points
        Wp(:,j) = [E(:,j); chi_m(5:7,j)-w_bar];

        % State Covariance Estimate
        Pm_p1 = Pm_p1 + 1/(2*n+1)*(Wp(:,j)*Wp(:,j).');
    end

%     P_iter = zeros(7,7);
%     Pm_p1 = zeros(7,7);
%     for j = 1:2*n+1
% 
%         %subtract quaternion components
%         chiC = EP2C(chi_m(1:4,j));
%         xmC = EP2C(xm_p1(1:4));
%         quat_C = chiC*xmC'; %transpose xmC?
%         quat_sub = C2EP(quat_C)/norm(C2EP(quat_C)); % Normalize error
% 
%         ang_sub = chi_m(5:7,j) - xm_p1(5:7); % Is this right?? ---
%         x_sub = [quat_sub ; ang_sub];
% 
%         % P_iter = w_c(:,j).*(chi_m(:,j) - xm_p1)*(chi_m(:,j) - xm_p1)' + Qk;
%         P_iter = w_c(:,j).*(x_sub)*(x_sub)' + Qk;
%         Pm_p1 = Pm_p1 + P_iter; 
%     end


%% 2. Measurement Update Step at time k+1 given observation y(k+1)
    %%%%%%%%%%%%%%%%%%%
    %part a - generate sigma pts
    %%%%%%%%%%%%%%%%%%%

    Sk_p1 = chol(Pm_p1,'lower');

    %repeat above steps for chi k+1
    chik0_p1 = xm_p1;

    for j = 1:n
        % W = sqrt(n+lambda)*Sk(j,:)';
        W = sqrt(2*n)*Sk_p1(j,:)';
        theta = norm(W(1:3)); % Principal Rotation of perturbations on attitude
        if theta < 1e-8
            q = [1 0 0 0].';
        else
            e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
            q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
            q = q/norm(q);
        end
        chik_p1(:,j) = [EP_Add(xm_p1(1:4),q); xm_p1(5:7)+W(4:6)]; % Combining the two
    end

    for j = (n+1):2*n
        % W = -sqrt(n+lambda)*Sk(j-n,:)';
        W = -sqrt(2*n)*Sk_p1(j-n,:)';
        theta = norm(W(1:3)); % Principal Rotation of perturbations on attitude
        if theta < 1e-8
            q = [1 0 0 0].';
        else
            e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
            q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
            q = q/norm(q);
        end
        chik_p1(:,j) = [EP_Add(xm_p1(1:4),q); xm_p1(5:7)+W(4:6)]; % Combining the two
    end

    chik_p1 = [chik0_p1 chik_p1]; %combine
    %%%%%%%%%%%%%%%%%%%
    %part b - propagate each chi through non-linear measurment function h
    %%%%%%%%%%%%%%%%%%%
    gam_p1 = zeros(3,2*n+1);

    for j = 1:2*n+1
        [gam_p1(:,j)] = h_ysim(chik_p1(:,j));
    end
    %%%%%%%%%%%%%%%%%%%
    %part c - get predicted measurement mean and measurement covar
    %%%%%%%%%%%%%%%%%%%   
    ym_p1(:,i) = 1/(2*n+1)*sum(gam_p1,2);

    Pyy_p1 = zeros(3);
    for j = 1:2*n+1
        P_iter_yy = 1/(2*n+1)*(gam_p1(:,j) - ym_p1(:,i))*(gam_p1(:,j) - ym_p1(:,i)).' + Rk;
        Pyy_p1 = Pyy_p1 + P_iter_yy; 
    end

    %%%%%%%%%%%%%%%%%%%
    %part d - get state measurement cross-covariance matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Pxy_p1 = zeros(6,3);
    for j = 1:2*n+1
    Pxy_p1_iter = 1/(2*n+1)*(Wp(:,j)*(gam_p1(:,j)-ym_p1(:,i)).');
    Pxy_p1 = Pxy_p1 + Pxy_p1_iter;
    end

    %%%%%%%%%%%%%%%%%%%
    %part e - Estimate kalman gain matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Kk_p1 = Pxy_p1*inv(Pyy_p1);

    %%%%%%%%%%%%%%%%%%%
    %part f - Perform kalman state and covariance update with observation yk+1 (nxp)
    %%%%%%%%%%%%%%%%%%%  
    update = Kk_p1*(y_t(:,i) - ym_p1(:,i));
    u_ang = norm(update(1:3));
    u_vec = update(1:3)/norm(update(1:3));
    u_quat = [cos(u_ang/2);u_vec*sin(u_ang/2)];
    u_quat = u_quat/norm(u_quat);
    
    xp_quat = EP_Add(xm_p1(1:4),u_quat);
    xp_quat = xp_quat/norm(xp_quat);

    xp_p1 = [xp_quat;xm_p1(5:7)+update(4:6)];
    Pp_p1 = Pm_p1 - Pxy_p1*inv(Pyy_p1)*Pxy_p1.';

    %Store variables
    x_ukf(:,i) = xp_p1;
    P_ukf(:,:,i) = Pp_p1;
    %Set for next iteration
    xp = xp_p1;
    Pp = Pp_p1;
    tlast = tvec(i);
end

NEES = 1;
NIS = 1; 

figure
plot(tvec,ym_p1(:,:))

end
