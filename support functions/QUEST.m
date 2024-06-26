function [b_av] = QUEST(b_vecs,w_vec)
%Note: b referese to quaternion vector, q refers to CRP vector

    B = zeros(3,3);
    for i = 1:length(w_vec)
        B = B + w_vec(i).*EP2C(b_vecs(:,i));
    end

    z = [(B(2,3) - B(3,2)) (B(3,1)-B(1,3)) (B(1,2)-B(2,1)) ]';
    S = B + B'; 
    K = [trace(B) z'
        z S-trace(B)*eye(3)    
        ];
    %Compute lambda using lambda characteristic polynomial:
    %Could also be computed less efficiently from 
    % K=[S-eye(3)*trB z; z' trB];
    [dummy, lam]=eigs(K, 1);

    % detS=det(S);s
    % adjS=detS*inv(S);
    % trB = trace(B); 
    % 
    % a=trB^2 - trace(adjS);  
    % b=trB^2 + z'*z;
    % c=detS + z'*S*z;
    % d=z'*S*S*z;
    % cnst=a*b + c*trB - d;
    % 
    % % lam=0.5*(sum(sum(qp.^2)) + sum(sum(pp.^2))); % 0.5*(trace(qp*qp') + trace(pp*pp'));
    % lam = sum(w_vec); %best guess for lambda
    % 
    % lamprev=0.0;
    % 
    % % Newton-Raphson
    % while abs((lam-lamprev)/lam) >= 1E-12
    %     lamprev=lam;
    %     lam=lam - (((lam^2 - (a+b))*lam - c)*lam + cnst)/(2*(2*lam^2 - (a+b))*lam - c);
    % end

    qbar = inv( (lam+trace(B))*eye(3) - S )*z;
    Bbar = (1/sqrt( (1+qbar'*qbar) )) * [1;qbar];
    quest_dcm = (Bbar(1)^2 - Bbar(2:4)'*Bbar(2:4))*eye(3) + 2*Bbar(2:4)*Bbar(2:4)' - 2*Bbar(1)*tilde(Bbar(2:4))
    b_av = C2EP(quest_dcm); 
end

