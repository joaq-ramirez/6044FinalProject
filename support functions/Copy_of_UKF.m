function [x_ukf,P_ukf,NEES,NIS] = UKF(x_t,y_t,tvec,c)
% load('orbitdeterm_finalproj_KFdata.mat')

xp = x_t(:,1); % initial value of x
Pp = 0.001*eye(6);  % Define covariance - 6 because only 6 DOF (quaternions are constrained to 3)
% Pp(1:3,1:3) = 0.01*eye(3);
% Rk = 0.001*eye(3); % Measurement noise - pre omega
Rk = 0.001*eye(3); % Measurement noise

Qk = zeros(6);
Qk(1:3,1:3) = 0.000001*eye(3); % Will need to adjust later Process Noise
Qk(4:6,4:6) = 0.000001*eye(3); % Will need to adjust later Process Noise

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
    i
    
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
            q = [1 0.5*W(1:3)']';
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
            q = [1 0.5*W(1:3)'].';
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
    omega_bar = 1/(2*n+1)*sum(chi_m(5:7,:),2);

    % Propogated State Estimate
    xm_p1 = [q_bar;omega_bar]; 

    % Estimating the State Covariance    
    Wp = zeros(6,length(chi_m(1,:)));
    Pm_p1 = zeros(6);

    for j = 1:length(chi_m(1,:))
        % Differences between average state and sigma points
        Wp(:,j) = [E(:,j); chi_m(5:7,j)-omega_bar];

        % State Covariance Estimate
        Pm_p1 = Pm_p1 + 1/(2*n+1)*(Wp(:,j)*Wp(:,j).');
    end
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
            q = [1 0.5*W(1:3)'].';
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
            q = [1 0.5*W(1:3)].';
        else
            e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
            q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
            q = q/norm(q);
        end
        chik_p1(:,j) = [EP_Add(xm_p1(1:4),q); xm_p1(5:7)+W(4:6)]; % Combining the two
    end

    if i == 345 
        zz = 1; 
    end

    chik_p1 = [chik0_p1 chik_p1]; %combine
    %%%%%%%%%%%%%%%%%%%
    %part b - propagate each chi through non-linear measurment function h
    %%%%%%%%%%%%%%%%%%%
    gam_p1 = zeros(3,2*n+1); % post omega

    for j = 1:2*n+1
        [gam_p1(:,j)] = h_ysim(chik_p1(:,j));
    end
    %%%%%%%%%%%%%%%%%%%
    %part c - get predicted measurement mean and measurement covar
    %%%%%%%%%%%%%%%%%%%   
    ym_p1(:,i) = 1/(2*n+1)*sum(gam_p1,2);
    ym_p1(1:3,i) = ym_p1(1:3,i)/ norm(ym_p1(1:3,i)); % normalize sun vector measure

    %Debug plot
    
%Debugging
    % if i == 500
    %     figure
    %     quiver3(zeros(1,13), zeros(1,13), zeros(1,13),gam_p1(1,:), gam_p1(2,:), gam_p1(3,:));
    %     hold on
    %     quiver3(0, 0,0,ym_p1(1,i), ym_p1(2,i), ym_p1(3,i),'r');
    %     quiver3(0, 0,0,y_t(1,i), y_t(2,i), y_t(3,i),'g');
    %     fprintf('debug')
    % 
    % end

    Pyy_p1 = zeros(3); % post omega
    for j = 1:2*n+1
        % meas_diff = (gam_p1(:,j) - ym_p1(:,i))/ norm(gam_p1(:,j) - ym_p1(:,i));
        % P_iter_yy = 1/(2*n+1)*(meas_diff)*(meas_diff');

        P_iter_yy = 1/(2*n+1)*(gam_p1(:,j) - ym_p1(:,i))*(gam_p1(:,j) - ym_p1(:,i)).';
        Pyy_p1 = Pyy_p1 + P_iter_yy; 
    end
    Pyy_p1 = Pyy_p1 + Rk; 
    %%%%%%%%%%%%%%%%%%%
    %part d - get state measurement cross-covariance matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Cxy_p1 = zeros(6,3); % post omega was 6,3
    for j = 1:2*n+1
        Pxy_p1_iter = 1/(2*n+1)*(Wp(:,j)*(gam_p1(:,j)-ym_p1(:,i))');
        Cxy_p1 = Cxy_p1 + Pxy_p1_iter;
    end

    %%%%%%%%%%%%%%%%%%%
    %part e - Estimate kalman gain matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Kk_p1 = Cxy_p1*inv(Pyy_p1);

    %%%%%%%%%%%%%%%%%%%
    %part f - Perform kalman state and covariance update with observation yk+1 (nxp)
    %%%%%%%%%%%%%%%%%%%  
    
    %Attempt 2
    % update = Kk_p1*(y_t(:,i) - ym_p1(:,i));
    % v_cross = cross(ym_p1(1:3,i), update(1:3));
    % w = v_cross/ norm(v_cross); 

    % norm_rot = w(1:3)/ norm(w(1:3)); 
    % theta_u = asin(min(max(norm(w), -1), 1)); 
    % theta_u = acos(dot(ym_p1(1:3,i), y_t(1:3,i))/ (norm(ym_p1(1:3,i))* norm(y_t(1:3,i))));
    % del_q =[cos(theta_u / 2); sin(theta_u / 2) * w];
    % q_update = quatmultiply(xm_p1(1:4)',del_q');
    % q_update = EP_Add(xm_p1(1:4),del_q);
    % q_update = q_update/ norm(q_update); 

    %Old attempt
    update = Kk_p1*(y_t(:,i) - ym_p1(:,i));
    u_ang = norm(update(1:3));
    u_vec = update(1:3)/norm(update(1:3));
    
    u_quat = [cos(u_ang/2);u_vec*sin(u_ang/2)];
    u_quat = u_quat/norm(u_quat);
    
    xp_quat = EP_Add(xm_p1(1:4),u_quat);
    % xp_quat = -(xm_p1(1:4).*u_quat);
    % xp_quat = xp_quat/norm(xp_quat);
    
    % xp_p1 = [q_update';xm_p1(5:7)+update(4:6)];
    xp_p1 = [xp_quat;xm_p1(5:7)+update(4:6)];
    % Pp_p1 = Pm_p1 - Cxy_p1*inv(Pyy_p1)*Cxy_p1';
    Pp_p1 = Pm_p1 - Kk_p1*Pyy_p1*Kk_p1'; 
    
    debug_update(:,i) = update; 
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
figure
plot(tvec,debug_update(1:3,:))
figure
plot(tvec,debug_update(4:6,:))

end
