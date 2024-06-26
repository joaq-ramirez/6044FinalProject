function [x_gsf,P_gsf] = GSF(x_t,y_t,tvec,c)

%% Gaussian Sum Filter
% Authors: Joaquin Ramirez, Alex Nelson

% Filter Component Initialization
xp = x_t(:,1); % Initial value of x
Pp = 0.01*eye(6);  % Define covariance - 6 because only 6 DOF (quaternions are constrained to 3)
Qk = 1e-3*eye(6); % Process Noise
wm = 1;
x_gsf = zeros(7,length(x_t));
P_gsf = zeros(6,6,length(x_t));
n = length(Pp(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Mixture Noise Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GMModel = GMM();

% Gaussian Sum Filter Loop
for k = 1:length(x_t)
    
    %% 1. dynamics prediction step from time step k->k+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run all Gaussian's of prior through UKF dynamics prediction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z = length(xp(1,:));
    for i = 1:z

        Sk = chol(Pp(:,:,i)+Qk,'lower');
        xm_p1 = zeros(7,z);
        Pm_p1 = zeros(6,6,z);
        wi = zeros(z);

        %%%%%%%%%%%%%%%%%%%
        %part A: generate sigma points from xp
        %%%%%%%%%%%%%%%%%%%
        chik0 = xp(:,i); %Initial chi_mean
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
            chik(:,j) = [EP_Add(xp(1:4,i),q); xp(5:7,i)+W(4:6)]; % Combining the two
        end
    
        for j = (n+1):2*n
            W = -sqrt(2*n)*Sk(j-n,:)';
            theta = norm(W(1:3)); % Principal Rotation of perturbations on attitude
            if theta < 1e-8
                q = [1 0 0 0].';
            else
                e = W(1:3)/norm(W(1:3)); % Principal Vector of perturbations on attitude
                q = [cos(theta/2);e*sin(theta/2)]; % Quaternion representation of perturbations 
                q = q/norm(q);
            end
            chik(:,j) = [EP_Add(xp(1:4,i),q); xp(5:7,i)+W(4:6)]; % Combining the two
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
    
        % Averaging the Quaternion portion of chi_m and recording error vectors
        % (for the Covariance Calcualtions)
        [q_bar, E] = Q_mean(chi_m(1:4,:));
    
        % Averaging the Angular Velocity Portion of chi_m
        w_bar = 1/(2*n+1)*sum(chi_m(5:7,:),2);
    
        % Propogated State Estimate
        xm_p1(:,i) = [q_bar;w_bar]; 
    
        % Estimating the State Covariance    
        Wp = zeros(6,length(chi_m(1,:)));
    
        for j = 1:length(chi_m(1,:))
            % Differences between average state and sigma points
            Wp(:,j) = [E(:,j); chi_m(5:7,j)-w_bar];
    
            % State Covariance Estimate
            Pm_p1(:,:,i) = Pm_p1(:,:,i) + 1/(2*n+1)*(Wp(:,j)*Wp(:,j).');
        end

        % Weight update (multiplied by one since Gaussian Process Noise)
        wi(i) = wm(i); 
    end

    %% 2. Measurement Update Step at time k+1 given observation y(k+1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Run all Gaussian's of prior through UKF Measurement Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z = length(xm_p1(1,:));
    v = length(GMModel.ComponentProportion);

    xp_p1 = zeros(7,z*v);
    Pp_p1 = zeros(6,6,z*v);
    wm = zeros(1,z*v);
    for s = 1:z
        %%%%%%%%%%%%%%%%%%%
        %part a - generate sigma pts
        %%%%%%%%%%%%%%%%%%%
        Sk_p1 = chol(Pm_p1(:,:,s),'lower');
    
        %repeat above steps for chi k+1
        chik0_p1 = xm_p1(:,s);
    
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
            chik_p1(:,j) = [EP_Add(xm_p1(1:4,s),q); xm_p1(5:7,s)+W(4:6)]; % Combining the two
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
            chik_p1(:,j) = [EP_Add(xm_p1(1:4,s),q); xm_p1(5:7,s)+W(4:6)]; % Combining the two
        end
    
        chik_p1 = [chik0_p1 chik_p1]; %combine
        %%%%%%%%%%%%%%%%%%%
        %part b - propagate each chi through non-linear measurment function h
        %%%%%%%%%%%%%%%%%%%
        gam_p1 = zeros(3,2*n+1);
    
        for j = 1:2*n+1
            [gam_p1(:,j)] = h_ysim(chik_p1(:,j));
        end

        for l = 1:v
            %%%%%%%%%%%%%%%%%%%
            %part c - get predicted measurement mean and measurement covar
            %%%%%%%%%%%%%%%%%%%   
            ym_p1 = 1/(2*n+1)*sum(gam_p1,2);
        
            Pyy_p1 = zeros(3);
            for j = 1:2*n+1
                Rk = GMModel.Sigma(:,:,l);
                P_iter_yy = 1/(2*n+1)*(gam_p1(:,j) - ym_p1)*(gam_p1(:,j) - ym_p1).' + Rk;
                Pyy_p1 = Pyy_p1 + P_iter_yy; 
            end
        
            %%%%%%%%%%%%%%%%%%%
            %part d - get state measurement cross-covariance matrix (nxp)
            %%%%%%%%%%%%%%%%%%%  
            Pxy_p1 = zeros(6,3);
            for j = 1:2*n+1
            Pxy_p1_iter = 1/(2*n+1)*(Wp(:,j)*(gam_p1(:,j)-ym_p1).');
            Pxy_p1 = Pxy_p1 + Pxy_p1_iter;
            end
        
            %%%%%%%%%%%%%%%%%%%
            %part e - Estimate kalman gain matrix (nxp)
            %%%%%%%%%%%%%%%%%%%  
            Kk_p1 = Pxy_p1*inv(Pyy_p1);
        
            %%%%%%%%%%%%%%%%%%%
            %part f - Perform kalman state and covariance update with observation yk+1 (nxp)
            %%%%%%%%%%%%%%%%%%%  
            mu = GMModel.mu(l,:).';
            update = Kk_p1*(y_t(:,k) - ym_p1-mu);
            u_ang = norm(update(1:3));
            u_vec = update(1:3)/norm(update(1:3));
            u_quat = [cos(u_ang/2);u_vec*sin(u_ang/2)];
            u_quat = u_quat/norm(u_quat);
            
            xp_quat = EP_Add(xm_p1(1:4,s),u_quat);
            xp_quat = xp_quat/norm(xp_quat);

            xp_p1(:,l*s) = [xp_quat;xm_p1(5:7,s)+update(4:6)];
            Pp_p1(:,:,l*s) = Pm_p1 - Pxy_p1*inv(Pyy_p1)*Pxy_p1.';
            wm(l*s) = wi(s)*GMModel.ComponentProportion(l);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ISSUE WITH ROUTINE (NEED TO DETERMINE HOW TO DO WEIGHTED QUATERNION
    % SUMMATION)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(wm)
    %Store variables
    x_gsf(:,k) = x_gsf(:,k)+wm(j)*xp_p1(:,j);
    P_gsf(:,:,k) = P_gsf(:,:,k)+wm(j)*(Pp_p1(:,:,j)+xp_p1(:,j)*xp_p1(:,j).');
    P_gsf(:,:,k) = P_gsf(:,:,k)-x_gsf(:,k)*x_gsf(:,k).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compress the Gaussian Mixtures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian Moment-Matching (Update in Future Versions of Filter)
xp = x_gsf(:,k);
Pp = P_gsf(:,:,k);
wm = 1;

end





