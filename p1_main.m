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

%Perfect measurement generation
% for i = 1:600
%     x(:,i+1) = non_linear_xsim(x(:,i),c);
% end
% N_Rs = [1; 0; 0];
% y_t = zeros(3,length(ss_meas));
% for i = 1:601
%     beta_BN = x(1:4,i);
%     BN = EP2C(beta_BN);
%     y_t(1:3,i) = BN*N_Rs;%+0.01*rand()*ones(3,1);
% end

%Mapping ss_meas to a sun vector for filter
% y_t = zeros(6,length(ss_meas));
y_t = zeros(3,length(ss_meas));
for i = 1:length(tvec)
    y_t(1:3,i) = pinv(ss_mapping)*ss_meas(i,:)';
    y_t(1:3,i) = y_t(1:3,i)/ norm(y_t(1:3,i));
    % y_t(4:6,i) = w_BN(i,:) + ones(1,3)*mvnrnd(0, 0.0000001); 
end


%% UKF Estimation

% [x_ukf,P_ukf,NEES,NIS] = UKF_6(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);
[x_ukf,P_ukf,NEES,NIS] = UKF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);

%Calculate without reflections
% ss_meas = load("datasets\Sun_Sensor_Data_Albedo.csv");
% y_t = zeros(3,length(ss_meas));
% for i = 1:length(tvec)
%     y_t(1:3,i) = pinv(ss_mapping)*ss_meas(i,:)';
%     y_t(1:3,i) = y_t(1:3,i)/ norm(y_t(1:3,i));
%     % y_t(4:6,i) = w_BN(i,:) + ones(1,3)*mvnrnd(0, 0.0000001); 
% end
% 
% [x_ukfr,P_ukf,NEES,NIS] = UKF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);


%% GSF Estimation
% [x_gsf,P_gsf] = GSF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);


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

% figure
% subplot(4,1,1)
% plot(tvec(3:601),x_ukfr(1,1:599),tvec(3:601),x_ukf(1,1:599),'LineWidth',2)
% title('UKF With Reflections v. Without Reflections')
% ylabel('\beta_0')
% subplot(4,1,2)
% plot(tvec(3:601),x_ukfr(2,1:599),tvec(3:601),x_ukf(2,1:599),'LineWidth',2)
% ylabel('\beta_1')
% subplot(4,1,3)
% plot(tvec(3:601),x_ukfr(3,1:599),tvec(3:601),x_ukf(3,1:599),'LineWidth',2)
% ylabel('\beta_2')
% subplot(4,1,4)
% plot(tvec(3:601),x_ukfr(4,1:599),tvec(3:601),x_ukf(4,1:599),'LineWidth',2)
% ylabel('\beta_3')
% legend(["Pred. without Reflections","Pred. with Reflections"])
% xlabel('Time (s)')

figure
subplot(4,1,1)
plot(tvec(3:601),x_t(1,3:601),'--',tvec(3:601),x_ukf(1,1:599),'LineWidth',2)
title('UKF Perfect Measure Quaternion State Estimate v. Truth')
ylabel('\beta_0')
subplot(4,1,2)
plot(tvec(3:601),x_t(2,3:601),'--',tvec(3:601),x_ukf(2,1:599),'LineWidth',2)
ylabel('\beta_1')
subplot(4,1,3)
plot(tvec(3:601),x_t(3,3:601),'--',tvec(3:601),x_ukf(3,1:599),'LineWidth',2)
ylabel('\beta_2')
subplot(4,1,4)
plot(tvec(3:601),x_t(4,3:601),'--',tvec(3:601),x_ukf(4,1:599),'LineWidth',2)
ylabel('\beta_3')
xlabel('Time (s)')


figure
subplot(3,1,1)
plot(tvec(3:601),x_t(5,3:601),'--',tvec(3:601),x_ukf(5,1:599),'LineWidth',2)
ylabel('\omega_1')
title('UKF Perfect Measure Angular Rate State Estimate v. Truth')
subplot(3,1,2)
plot(tvec(3:601),x_t(6,3:601),'--',tvec(3:601),x_ukf(6,1:599),'LineWidth',2)
ylabel('\omega_2')
subplot(3,1,3)
plot(tvec(3:601),x_t(7,3:601),'--',tvec(3:601),x_ukf(7,1:599),'LineWidth',2)
ylabel('\omega_3')
xlabel('Time (s)')

