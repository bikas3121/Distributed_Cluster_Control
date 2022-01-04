function [valP, valF] = solverLMI(nP, sigma, A_k_cl, Q_k)
%Solves the LMI and returns the value of P and the LMI
% if the eigenvlaues of the P are all positives and the eigenvalues of F
% are all negatives then the LMI condition in the paper is satisfied.

%% Inputs:
% nP - dimenstion of the matrix variable or the decision variable.
% sigma - user input for the bounds. 
% The LMI has the form:
% F : A^TP + PA + Q_k < =0
% P < sigma*I_np
% P >0 - P is positive definite

%% Outputs:
% valP - matrix P 
% valF - LMI


epsilon=10^-3; %tolerance for the strict inequality
 

%sigmaM = sigma*eye(nP);
P = sdpvar(nP,nP);
F = [P >= epsilon*eye(nP), P <= sigma*eye(nP),  A_k_cl'*P + P*A_k_cl + Q_k  <= -epsilon*eye(nP)];

%F = [P >= epsilon*eye(nP), P <= sigma*eye(nP),  A_k_cl'*P + P*A_k_cl  <= -epsilon*eye(nP)];

options = sdpsettings('solver','sedumi');
sol = optimize(F, -sigma, options);
feas = sol.problem;
if feas == 1 
    disp('************************************')
    disp('Not feasible')
    disp('************************************')
else
    disp('************************************')
    disp('Feasible')
    disp('************************************')
end
valP = value(P);

valF = A_k_cl'*valP + valP*A_k_cl + Q_k;
end