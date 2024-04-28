%% Reflections
% Adds in Reflections Randomly to the Sun Sensor Measurements from Basilisk

clear
close all
clc

% Pull in CSS data set
dat = load("..\datasets\Sun_Sensor_Data_Albedo.csv");

% Set Random Number Generator
rng(22)

% Probability of a reflection
p = 0.01;

for j = 1:20
    for i = 1:length(dat)
        if rand() < p
            dat(i,j) = dat(i,j) + (3-dat(i,j))*rand();
        end
    end
end 

csvwrite("..\datasets\Sun_Sensor_Data_Albedo_Reflections.csv",dat);

dat = load("..\datasets\Sun_Sensor_Data_Albedo_Reflections.csv");

% Plot Results to check outputs
t = 0:1:600;
figure(1)
plot(t,dat(:,1),...
    t,dat(:,2),...
    t,dat(:,3),...
    t,dat(:,4),...
    t,dat(:,5),...
    t,dat(:,6),...
    t,dat(:,7),...
    t,dat(:,8),...
    t,dat(:,9),...
    t,dat(:,10),...
    t,dat(:,11),...
    t,dat(:,12),...
    t,dat(:,13),...
    t,dat(:,14),...
    t,dat(:,15),...
    t,dat(:,16),...
    t,dat(:,17),...
    t,dat(:,18),...
    t,dat(:,19),...
    t,dat(:,20));

figure(2)
subplot(5,1,1)
plot(t,dat(:,1))
title("Sun Sensor Readings with Reflections")
ylabel("CSS1 Voltage")
subplot(5,1,2)
plot(t,dat(:,2))
ylabel("CSS2 Voltage")
subplot(5,1,3)
plot(t,dat(:,3))
ylabel("CSS3 Voltage")
subplot(5,1,4)
plot(t,dat(:,4))
ylabel("CSS4 Voltage")
subplot(5,1,5)
plot(t,dat(:,5))
ylabel("CSS5 Voltage")
xlabel("time [sec]")
