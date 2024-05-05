function [a_tilde] = tilde(a)
%tilde(a) Returns the tilde matrix of the vector 'a' (3 dim)

% Throw error if vector has more than 3 elements
if numel(a) > 3
    error('Your vector is too long')
end

% Decompose vector into elements
a1 = a(1); a2 = a(2); a3 = a(3);

% Construct tilde matrix
a_tilde = [ 0 -a3 a2; a3 0 -a1; -a2 a1 0];
end