% non_linear_x_sim_test

clear
close all
clc

t = 1:1:100;
x = zeros(7,100);
c.I = [0.012 0 0;0 0.01 0;0 0 0.006]; 
x(:,1) = [0.754385964912281;
      0.175438596491228;
      0.350877192982456;
      -0.526315789473684;
      0.0174532925199433;
      0.0349065850398866;
      -0.0174532925199433];

for i = 1:99
    x(:,i+1) = non_linear_xsim(x(:,i),c);
end

plot(t,x(1,:),t,x(2,:),t,x(3,:),t,x(4,:),t,x(5,:),t,x(6,:),t,x(7,:))