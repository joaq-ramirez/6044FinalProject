function [GMModel] = GMM()
    %% Gaussian Mixture Model of Sensor Noise
    
    % 
    % This function generates a mixed gaussian estimation of the sensor noise
    % based on data loaded in from the Basilisk Simulation.
    %
    
    % Generating Measurement Noise Data from Basilisk Sim Results
    dat = load(".\datasets\GMM_Sun_Sensor_Data_No_Noise.csv");
    datn = load(".\datasets\GMM_Sun_Sensor_Data_Noise_Plus.csv");
    
    % Converting Sun Sensor Measurements to Sun Estimation Vectors
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
    
    y = zeros(3,length(dat(:,1)));
    yn = zeros(3,length(datn(:,1)));
    for i = 1:length(datn(:,1))
        y(:,i) = pinv(ss_mapping)*dat(i,:)';
        y(:,i) = y(:,i)/norm(y(:,i));
        yn(:,i) = pinv(ss_mapping)*datn(i,:)';
        yn(:,i) = yn(:,i)/norm(yn(:,i));
    end
    
    noise = yn.'-y.';
    noise = noise(3:length(noise),:); % Truncated becuase index 1 is NAN
    
    % Creating a Gaussian Mixture Model of the Sensor Noise
    options = statset('MaxIter',1000);
    GMModel = fitgmdist(noise,5,'options',options);
end
 
