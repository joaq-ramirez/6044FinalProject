%% UKF Function
%Authors: Joaquin Ramirez, Alex Nelson

%% Housekeeping
clear all; clc; close all; 
rng(100)
%% Load in test data set to validate UKF
load('orbitdeterm_finalproj_KFdata.mat')

rE = 6378;
mu = 4*10^5;
omega_e = 2*pi/86400;

perturb_x0 = [0,0.075,0,-0.021]';
x0 = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';

B = [0 0;
    1 0;
    0 0;
    0 1];
Gamma = [0 0;
        1 0;
        0 0;
        0 1];
delT = 10;
x = x0+perturb_x0;
P = 2*pi*sqrt(6678^3/mu);

%data vars: Qtrue, Rtrue, tvec, ydata
%Important Notes:
%No input u
dt = delT;
Qk = Qtrue; % Qtrue from provided data set


Omk = dt*[0 0 ; 1 0; 0 0 ; 0 1];
Sw = chol(Qk,'lower');
Sv = chol(Rtrue,'lower'); % Rtrue from provided data set

% Monte Carlo Sim

N = 10;
n = 2;

alpha = 0.05;

for m = 1:N

    x_k = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';

%Noisy Data sim of x
for i = 1:length(tvec)

    %Simulate truth
    w_k = mvnrnd(zeros(2,1),Qk); % of Qtrue
    t_k = tvec(i);

    [~,x_nom] = ode45(@(t,x) odeLKF(t_k,x,[0;0]),[0 dt],x_k);

    x_nom = x_nom'; 
    
    x_k = x_nom(:,end);

    xnomp(:,i) = x_nom(:,end);

    w = Sw*randn(2,1);
    x_t(:,i) = xnomp(:,i) + Omk*w; % Truth solution
    %Simulate measurement truth
    for j = 1:12
        v = Sv*randn(3,1);
        y_t(:,i,j) = dynamicsCalc(x_t(:,i),j,t_k)+v;
        if any(~isnan(y_t(:,i,j)))
            station(:,i,1) = j;            
        end
    end
end

[x_ukf,P_ukf,NEES,NIS] = UKF(x_t,y_t,station,tvec);

xdiff = x_t - x_ukf;

end

NEESmean = mean(NEES);
NISmean = mean(NIS);

function [x_ukf,P_ukf,NEES,NIS] = UKF(x_t,y_t,station,tvec)
load('orbitdeterm_finalproj_KFdata.mat')
mu = 398600;

% y_d = y_t()
% if any(any(isnan(y_d))) || any(any(isempty(y_d)))
% %             yk = ystat(:,:,station(i));
% %             yk = yk(:,i);
%             stationCur = station(i);
%         else
%             stationCur = y_d(4);
%             yk = y_d(1:3,1);
%         end

dt = 10;
Omk = dt*[0 0 ; 1 0; 0 0 ; 0 1];
perturb_x0 = [0,0.075,0,-0.021];

P0 = 1e-3*eye(4);
P_p = P0; %What is P0? set it smalll..
dx_p = perturb_x0';
Qk = Qtrue; 
% Qk = [0.02 0; 0 0.02];
Rp1 = Rtrue; 

xp = [6678, 0, 0, 6678 * sqrt(mu/(6678^3))]';
Pp = 100*eye(4,4);
y_k = zeros(3,1);
Qk = Qtrue;
Rk = Rtrue;

%preallocate for UKF loop
n = length(x_t(:,1)); 
kappa = 0;
beta = 2; 
alpha = 1e10;% estimate
lambda = alpha^2 * (n+kappa)-n; 
tlast = 0;

for i = 1:length(tvec)
    %reset each loop
%     Pp = 100*eye(4,4);
    Sk = chol(Pp,'lower');    
    Qk = eye(4); %Will need to adjust later
%     Pm_p1 = zeros(4,4);
    Pyy_p1 = zeros(3,3);
    Pxy_p1 = zeros(4,3);

    
%1. dynamics prediction step from time step k->k+1
    %%%%%%%%%%%%%%%%%%%
    %part a
    %%%%%%%%%%%%%%%%%%%
    chik0 = xp;
    for j = 1:n
        chik(:,j) = xp + (sqrt(n+lambda)*Sk(j,:)');
    end
    for j = (n+1):2*n
        chik(:,j) = xp - (sqrt(n+lambda)*Sk(j-n,:)');
    end
    chi_comb = [chik0, chik];
    inp_ic = chi_comb; 
    %%%%%%%%%%%%%%%%%%%
    %part b
    %%%%%%%%%%%%%%%%%%%
    for j = 1:2*n+1
        func = @(t,inp)[inp(2); -mu*inp(1)/(sqrt((inp(1)^2+inp(3)^2))^3);
            inp(4); -mu*inp(3)/(sqrt((inp(1)^2+inp(3)^2))^3)];
        tspan = [0 10]; %%%%%%%%%%%
        [~, xout] = ode45(func,tspan,inp_ic(:,j));
        chi_m(:,j) = xout(end,:)';
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
    xm_p1 = sum(w_m.*chi_m,2); % sum along dim 2
    
    P_iter = zeros(4,4);
    Pm_p1 = zeros(4,4);
    for j = 1:2*n+1
        P_iter = w_c(:,j).*(chi_m(:,j) - xm_p1)*(chi_m(:,j) - xm_p1)' + Qk;%THIS IS THE ERROR
        Pm_p1 = Pm_p1 + P_iter; 
    end
%2. Measurement Update Step at time k+1 given observation y(k+1)
    %%%%%%%%%%%%%%%%%%%
    %part a - generate sigma pts
    %%%%%%%%%%%%%%%%%%%
    Sk_p1 = chol(Pm_p1,'lower');
%     Sk_p1 = chol(eye(4),'lower');
    %repeat above steps for chi k+1
    chik0_p1 = xm_p1;
    for j = 1:n
        chik_p1(:,j) = xm_p1 + (sqrt(n+lambda)*Sk_p1(j,:)');
    end
    for j = (n+1):2*n
        chik_p1(:,j) = xp - (sqrt(n+lambda)*Sk_p1(j-n,:)');
    end
    chik_p1 =[chik0_p1 chik_p1]; %combine
    %%%%%%%%%%%%%%%%%%%
    %part b - generate sigma pts
    %%%%%%%%%%%%%%%%%%%    
    for j = 1:2*n+1
        [~,gam_p1(:,j),~,~] = linearizeH(i,chik_p1(:,j),tvec,station(:,i));
    end
    %%%%%%%%%%%%%%%%%%%
    %part c - get predicted measurement mean and measurement covar
    %%%%%%%%%%%%%%%%%%%   
    ym_p1 = sum(w_m.*gam_p1,2);
    P_iter_yy = zeros(3,3);
    for j = 1:2*n+1
        P_iter_yy = w_c(:,j).*(gam_p1(:,j) - ym_p1)*(gam_p1(:,j) - ym_p1)' + Rk;
        Pyy_p1 = Pyy_p1 + P_iter_yy; 
    end
    %%%%%%%%%%%%%%%%%%%
    %part d - get state measurement cross-covariance matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    for j = 1:2*n+1
    Pxy_p1_iter = w_c(:,j).*(chik_p1(:,j)- xm_p1)*(gam_p1(:,j)-ym_p1)'; %%
    Pxy_p1 = Pxy_p1 + Pxy_p1_iter;
    end
    %%%%%%%%%%%%%%%%%%%
    %part e - Estimate kalman gain matrix (nxp)
    %%%%%%%%%%%%%%%%%%%  
    Kk_p1 = Pxy_p1*inv(Pyy_p1);
    %%%%%%%%%%%%%%%%%%%
    %part f - Perform kalman state and covariance update with observation yk+1 (nxp)
    %%%%%%%%%%%%%%%%%%%    
    xp_p1 = xm_p1 + Kk_p1*(y_t(:,i,station(i))- ym_p1);
    Pp_p1 = Pm_p1 - Pxy_p1*inv(Pyy_p1)*Pxy_p1';
    
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

end


