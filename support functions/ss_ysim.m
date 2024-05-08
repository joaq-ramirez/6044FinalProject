function [gamma_m] = ss_ysim(chik_p1)
%non_linear_ysim Map x state to y measuremnts of the sun vector

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

gamma_m = ss_mapping*chik_p1(1:3); 

% for i = 1:length(ss_mapping)
%     gamma_m(:,i) =dot(ss_mapping(i,:),chik_p1(1:3))/(norm(ss_mapping(i,:))*norm(chik_p1(1:3)));
% end
% gamma_m(gamma_m>1)= nan; 

gamma_m(gamma_m<0) = 0; 
gamma_m = 3*gamma_m; 

end