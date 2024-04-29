function [q_bar, E] = Q_mean(Q,q,tol)
    
% Q_mean(Q,q,tol)
% 
%   q_bar = Q_mean(Q,q,tol) returns the mean quaternion given a set of
%   quaternions contained in Matrix Q. The fucntion is intialized with
%   quaternion guess, q, and then the error vector between the quaternions
%   in Q and the guess q is calculated. The sum of these error vectors are
%   summed to give a guess of q_bar. If the resulting error vector has a
%   magnitude over the tolerance, q = q_bar and the calculations are
%   repeated.
%

e = [0;0;0];
E = zeros(3,length(Q(1,:)));
error = 1;

while error > tol

    for i = 1:length(Q(2,:))
        qi = Q(:,i); % quaternion associated with index i
    
        q_inv = [q(1); -q(2); -q(3); -q(4)]; % inverse rotation of quaternion estimate
        
        qe = EP_Add(q_inv,qi); % Error Rotation between q attitude and qi (quaternion subtraction)
        
        % check that the error rotation is correct
        qi_test = EP_Add(q,qe);
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
    
        e = e+1/(2*6)*ei;

        E(:,i) = ei;
    end
    
    a = norm(e);
    v = e/norm(e);
    err = [cos(a/2);v*sin(a/2)];
    
    q = EP_Add(q,err);
    
    error = norm(e);
end

q_bar = q;
end