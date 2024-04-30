%% UKF Function
%Authors: Joaquin Ramirez, Alex Nelson

%% Housekeeping
clear all; clc; close all; 
rng(100)
addpath 'C:\Users\alexn\Documents\Graduate School\Spring 2024\ASEN 6044\6044FinalProject\support functions'

%% Data load
%Load quaternion and w data
B_BN = load("datasets\EP_Attitude_BN.csv");
w_BN = load("datasets\Angular_Velocity_BN_B.csv");
x_t = [B_BN, w_BN]'; 
tvec = 0:1:10*60;

%load sun sensor data
ss_meas = load("datasets\Sun_Sensor_Data_Albedo.csv");
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
y_t = zeros(3,length(ss_meas));
for i = 1:length(tvec)
    y_t(:,i) = inv(ss_mapping.'*ss_mapping)*ss_mapping.'*ss_meas(i,:)';
end

plot(tvec,y_t)


%% Declare constants
c.I = [0.012 0 0;0 0.01 0;0 0 0.006]; 

% EP Attitude [0.7543859649122806, 0.087719298245614, 0.1754385964912281, -0.2631578947368421] (B/N)
% Inertial Position [4863577.031787383, 4863577.031787383, 0.0] (m in inertial frame)
% Inertial Velocity [5382.927016299871, -5382.927016299871, 0.0] (m/s in inertial frame)
% Angular Velocity [0.017453292519943295, 0.03490658503988659, -0.017453292519943295] (rad/s of B/N in B frame)

[x_ukf,P_ukf,NEES,NIS] = UKF(x_t(:,3:length(x_t)),y_t(:,3:length(y_t)),tvec(3:length(tvec)),c);

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

figure(1)
subplot(4,1,1)
plot(tvec(3:601),x_t(1,3:601),tvec(3:601),x_ukf(1,1:599))
subplot(4,1,2)
plot(tvec(3:601),x_t(2,3:601),tvec(3:601),x_ukf(2,1:599))
subplot(4,1,3)
plot(tvec(3:601),x_t(3,3:601),tvec(3:601),x_ukf(3,1:599))
subplot(4,1,4)
plot(tvec(3:601),x_t(4,3:601),tvec(3:601),x_ukf(4,1:599))

figure(2)
subplot(3,1,1)
plot(tvec(3:601),x_t(5,3:601),tvec(3:601),x_ukf(5,1:599))
subplot(3,1,2)
plot(tvec(3:601),x_t(6,3:601),tvec(3:601),x_ukf(6,1:599))
subplot(3,1,3)
plot(tvec(3:601),x_t(7,3:601),tvec(3:601),x_ukf(7,1:599))
