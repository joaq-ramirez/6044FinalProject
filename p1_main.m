%% UKF Function
%Authors: Joaquin Ramirez, Alex Nelson

%% Housekeeping
clear all; clc; close all; 
rng(100)
addpath '.\support functions'

%% Data load
%Load quaternion and w data
B_BN = load("datasets\EP_Attitude_BN.csv");
w_BN = load("datasets\Angular_Velocity_BN_B.csv");
x_t = [B_BN, w_BN]'; 
tvec = 0:1:10*60;

%load sun sensor data
% ss_meas = load("datasets\Sun_Sensor_Data_Albedo.csv");
ss_meas = load("datasets\Sun_Sensor_Data_Albedo_Reflections.csv");
% ss_meas = load("datasets\Sun_Sensor_Data_No_Albedo.csv");

ss_mapping = [-0.866,-0.500,0.000;
              -0.866,0.000,0.500;
              -0.866,0.500,0.000;
              -0.866,0.000,-0.500;
              0.000,0.500,0.866;
              -0.500,0.000,0.866;
              0.000,-0.500,0.866;
              0.500,0.000,0.866;
              0.000,0.500,-0.866;
              0.500,0.000,-0.866;
              0.000,-0.500,-0.866;
              -0.500,0.000,-0.866;
              -0.500,0.866,0.000;
              0.000,0.866,0.500;
              0.500,0.866,0.000;
              0.000,0.866,-0.500;
              0.000,-0.866,-0.500;
              0.500,-0.866,0.000;
              0.000,-0.866,0.500;
              -0.500,-0.866,0.000];



%% Declare constants
c.I = [0.012 0 0;0 0.01 0;0 0 0.006]; 

x = zeros(7,100);
x(:,1) = [0.754385964912281;
      0.175438596491228;
      0.350877192982456;
      -0.526315789473684;
      0.0174532925199433;
      0.0349065850398866;
      -0.0174532925199433];

for i = 1:600
    x(:,i+1) = non_linear_xsim(x(:,i),c);
end

N_Rs = [1; 0; 0];
y_t = zeros(3,length(ss_meas));
for i = 1:601
    beta_BN = x(1:4,i);
    BN = EP2C(beta_BN);
    y_t(1:3,i) = BN*N_Rs;%+0.01*rand()*ones(3,1);
end

for i = 1:length(ss_meas)
    noisy_y(:,i) = pinv(ss_mapping)*ss_meas(i,:)';
end

% %Add in angular rates
% % y_t = zeros(6,length(ss_meas));
% for i = 1:length(tvec)
%     % y_t(1:3,i) = pinv(ss_mapping)*ss_meas(i,:)';
%     % y_t(1:3,i) = y_t(1:3,i)/ norm(y_t(1:3,i));
%     y_t(4:6,i) = w_BN(i,:) ;%+ ones(1,3)*mvnrnd(0, 0.0000001); 
% end
% 
% figure
% plot(tvec,y_t(1:3,:))
% figure
% plot(tvec,y_t(4:6,:))

% plot(tvec,y_t)

% EP Attitude [0.7543859649122806, 0.087719298245614, 0.1754385964912281, -0.2631578947368421] (B/N)
% Inertial Position [4863577.031787383, 4863577.031787383, 0.0] (m in inertial frame)
% Inertial Velocity [5382.927016299871, -5382.927016299871, 0.0] (m/s in inertial frame)
% Angular Velocity [0.017453292519943295, 0.03490658503988659, -0.017453292519943295] (rad/s of B/N in B frame)

%% UKF Estimation

% [x_ukf,P_ukf,NEES,NIS] = UKF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);

%% GSF Estimation
% [x_gsf,P_gsf] = GSF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);

%% Sun Sensor UKF
ys_t = ss_meas(3:end,:); 
xs_t = [y_t; w_BN'] ;
[x_ukf,P_ukf,NEES,NIS] = UKFSunSensor(xs_t,ys_t,tvec(3:length(ys_t)),c); % y_t is our new state

%% Inital plotting of sim data

% figure
% ylabl = ["\Beta_1", "\Beta_2", "\Beta_3", "\omega_1", "\omega_2", "\omega_3"]; 
% sgtitle('Mission Data Sim')
% pos_v = [1 3 5 2 4 6];
% for i = 1:6
%     subplot(3,2,pos_v(i))
%     plot(tvec,x_t(i+1,:),'linewidth',2)
%     xlabel("Time [s]")
%     ylabel(ylabl(i))
% end

figure
subplot(3,1,1)
title('UKF Sun Sensor Pred. w/Noisy/Reflection Sensor Data')
plot(tvec(3:length(ys_t)),xs_t(1,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(1,1:length(ys_t)-2),tvec(3:601),noisy_y(1,1:599))
ylabel('Sun Bector {}^Bx_1')
subplot(3,1,2)
plot(tvec(3:length(ys_t)),xs_t(2,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(2,1:length(ys_t)-2),tvec(3:601),noisy_y(2,1:599))
ylabel('Sun Bector {}^Bx_2')
subplot(3,1,3)
plot(tvec(3:length(ys_t)),xs_t(3,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(3,1:length(ys_t)-2),tvec(3:601),noisy_y(3,1:599))
ylabel('Sun Bector {}^Bx_3')
legend(["Truth","UKF Prediction", "Noisy Data"])
xlabel('Time (s)')
for i = 1:length(ys_t)-2
   ang(i) =  acos(dot(xs_t(:,i+2), x_ukf(:,i)));
end

figure
plot(tvec(3:length(ys_t)),ang)


figure
subplot(3,1,1)
title('UKF Angular Rate Pred. w/Noisy/Reflection Sensor Data')
plot(tvec(3:length(ys_t)),xs_t(4,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(4,1:length(ys_t)-2))
ylabel('{}^B\omega_1')
subplot(3,1,2)
plot(tvec(3:length(ys_t)),xs_t(5,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(5,1:length(ys_t)-2))
ylabel('{}^B\omega_2')
subplot(3,1,3)
plot(tvec(3:length(ys_t)),xs_t(6,3:length(ys_t)),tvec(3:length(ys_t)),x_ukf(6,1:length(ys_t)-2))
ylabel('{}^B\omega_3')
legend(["Truth","UKF Prediction", "Noisy Data"])
xlabel('Time (s)')
