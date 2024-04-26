%% UKF Function
%Authors: Joaquin Ramirez, Alex Nelson

%% Housekeeping
clear all; clc; close all; 
rng(100)

%% Data load
%Load quaternion and w data
B_BN = load("datasets\EP_Attitude_BN.csv");
w_BN = load("datasets\Angular_Velocity_BN_B.csv");
x_t = [B_BN, w_BN]'; 
tvec = 0:1:10*60;

%load sun sensor data
y_meas = load("datasets\Sun_Sensor_Data_Albedo.csv");
ss_mapping = ones(3,20); % Placeholder
for i = 1:length(tvec)
    y_t(:,i) = ss_mapping*y_meas(i,:)';
end



%% Declare constants
c.I = eye(3).*[10 50 20]; 

% EP Attitude [0.7543859649122806, 0.087719298245614, 0.1754385964912281, -0.2631578947368421] (B/N)
% Inertial Position [4863577.031787383, 4863577.031787383, 0.0] (m in inertial frame)
% Inertial Velocity [5382.927016299871, -5382.927016299871, 0.0] (m/s in inertial frame)
% Angular Velocity [0.017453292519943295, 0.03490658503988659, -0.017453292519943295] (rad/s of B/N in B frame)

[x_ukf,P_ukf,NEES,NIS] = UKF(x_t,y_t,tvec,c);

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