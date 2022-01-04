function [Kext_opt] = suboptimal_gain_generation_algorithm(nCk,a, b2, Lext, init, Kext, Rk, sigma)

% The function [Kext_opt] = suboptimal_gain_generation_algorithm(nCk,a, b2,
% Lext, init, Kext, Rk, sigma) generates the suboptimal gain that satisfies
% the LMI conditions in the equation 43 in the paper. Starting from the
% initial gains that synchronize the system it calculates the suboptimal
% gain profile that satisfies the satisfaction equillibrium. 

% INPUTS:
% nCk - vectors of the agents cardinality in the network. 
% a - agent's state matrix. 
% b2 - agent's input matrix. 
% Lext - External Laplcian matrix of the network without internal
% connections. 
% init - initial conditions. 
% Kext - initial gain profile that synchronize the system. 
% Rk - weighing matrix. 
% sigma - given threshold. 

% OUTPUTS:
% Kext_opt - The suboptimal gain that satisfies the satisfaction
% equillibrium condition in the equation 43 in the paper. 




 %%
 nx = 2;
 [val_count] = initial_gain_check(nCk, nx, a, b2, Lext, init, Kext, Rk, sigma);
 if val_count  == length(nCk)
       disp('********************************************')
    disp('Good Initial Gain Profile ')
    disp('********************************************')
 else
       disp('********************************************')
    disp('Bad Initial Gain Profile')
    disp('********************************************')
 end
 
 pause(5)
% dbstop in main2_alg.m at 96 if (val_count ~= length(nCk))   

%% Optimal Gain generation

% Cluster 1;
C1 = 1;
[Kopt1, ~,~,~] = sub_opt_gain (nCk, nx, Kext, C1, a, b2, Lext, init, Rk, sigma);
C2 = 2;
[Kopt2, ~,~,~] = sub_opt_gain (nCk, nx, Kopt1, C2, a, b2, Lext, init, Rk, sigma);
C3 = 3;
[Kopt3, ~,~,~] = sub_opt_gain (nCk, nx, Kopt2, C3, a, b2, Lext, init, Rk, sigma);
C4 = 4;
[Kopt4, ~,~,~] = sub_opt_gain (nCk, nx, Kopt3, C4, a, b2, Lext, init, Rk, sigma);

Kext_opt = Kopt4;

end