function [Kint_opt, Kext] = initial_gain(nCk, a, b2, Lextbar, Rk)

% The function [Kint_opt, Kext] = initial_gain(nCk, a, b2, Lextbar, Rk)
% generatest internal and the external gains for the given network
% dynamics. The internal gain is calculated using the relation in Lemma 2
% in the design procedure in the paper. The initial external gain profile
% synchronizing the system can be choosen as random. We use the algorithm
% provided in the article "Guaranteed cost control design for
% synchronization in netowrks of linearly singularly perturbed system,
% J.B.Rejeb, IC Morarescu and J. Daafouz". 

% Inputs:
% nCk  - vector containing the number of agents in the network
% a - state matrix of agent
% b2 - input matrix 
% Lextbar - the weighted Laplacian matrix. See equation Lextbar below
% equation 33 in the paper. 
% Rk - is the symmetric weighing matrix. 

% Outputs:
% Kint_opt - The optimal internal gain obtained using LQR 
% Kext - Initial external gain profile synchronizing the network . 

%% Test Parameters
% nCk = [5 10 15 20];
% 
% theta = 30; 
% a =1*[cos(theta),-sin(theta);sin(theta),cos(theta)];
% b2=[1;1];  
% Calculate the Gain -LQR

%% Internal Gain 
nx = 2;
m = length(nCk);
% Rk = [];
% for i = 1:m
% %     Rk =[Rk,  eye(1)];
%     Rk = [Rk, eye(1)*(nCk(i)^(-1))];
% end

Kint = [];
for i = 1:m
    Qint=nCk(i)*eye(nx);
    Rint=nCk(i)^2*Rk(i);
    [Kint1,Pint,e] = lqr(a,nCk(i)*b2,Qint,Rint);
    Kint = [Kint; Kint1];
end
Kint_opt = 1*Kint;

%% External Gain
% External Gain for all clusters using Jihenes algorithm
Kext = [];
for i = 1:m
    R = Rk(i);
%     R = 1;
    [P_eps,Kextj,G] = Gain_Ji_Jihene(a, b2, Lextbar,R);
    Kext = [Kext; 1*Kextj];
end
  
end