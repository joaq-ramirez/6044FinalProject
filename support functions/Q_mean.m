function [q_bar, E] = Q_mean(Q)
    
% Q_mean(Q,q,tol)
% 
%   [q_bar, E] = Q_mean(Q) returns the mean quaternion (q_bar) given a set
%   of quaternions contained in matrix Q. It also returns a matrix (E) of 
%   the error vectors of each quaternion component of Qi compared to the 
%   mean rotation.
%

E = zeros(3,length(Q(1,:)));

quat = quaternion(Q.');

q_bar = compact(meanrot(quat)).';

for i = 1:length(Q(1,:))
    qi = Q(:,i); % quaternion associated with index i

    q_inv = [q_bar(1); -q_bar(2); -q_bar(3); -q_bar(4)]; % inverse rotation of quaternion estimate
    
    qe = EP_Add(q_inv,qi); % Error Rotation between q attitude and qi (quaternion subtraction)
    
    % check that the error rotation is correct
    qi_test = EP_Add(q_bar,qe);
    if norm(qi_test-qi) > 1e-4
        fprintf("error computing quaternion inverse:"+string(i)+"\n")
    end

    theta = 2*acos(qe(1));
    if theta < 1e-5
        ei = [0;0;0];
    else
        e1 = qe(2)/sin(theta/2);
        e2 = qe(3)/sin(theta/2);
        e3 = qe(4)/sin(theta/2);
        ei = [e1;e2;e3];
    end

    E(:,i) = ei;
end
end