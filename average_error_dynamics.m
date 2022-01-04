function  [t_avg, Y_k, A_k_k,  J_avg_inst] = average_error_dynamics  (nCk, a, b2, Ck, Kext, Lext,  Rk, init, tspan, delta_t)

% The function  [t_avg, Y_k, A_k_k,  J_avg_inst] = average_error_dynamics 
% (nCk, a, b2, Ck, Kext, Lext,  Rk, init, tspan, delta_t) calcutes the
% average error dynamics in the equation 38 in the paper. 


% INPUTS:
% nCk - vectors containing the number of agents in the cluster.
% Ck - represents the corresponding cluster.
% a - internal dynamics : state matrix of an agent
% b2 - internal dynamics : input matrix of an agent
% Lext - External Laplacian of the network 
% Kint - The optimal  internal gain obtained from LQR
% Kext - External gain that satisfies the LMI (eq 48 in the paper)
% Rk - Weighing matrix in the original cost 
% init & tspan - Initial values and the time span for the solving the odes

% OUTPUTS:
% t_avg - time overwhich the dynamics is solved
% Y_k - state value of the average error dynamics
% A_k_k -  the state matrix of the average error dynamics in the equation
% 38
% J_avg_inst - Instantaneous average cost . 


%%
nx =2;
m = length(nCk);
[~,  H, W] = weighted_laplacian(nCk, Lext);

%%  Initial Value Generation
% Geneates the initial error condition for the given Ck                                         
initavg  =kron(W, eye(2))*kron(H,eye(2))'*init;  % average of the initial values (init), which is the initial conditions for the bar_y variable
% initavg
% Initial Values for the average errro dynamics=
if Ck == 1
    init_Ck = initavg(1:2,:);
    init_Ck_minus = initavg;
    init_Ck_minus(1:2,:) = [];
else
    st  = Ck*nx-1;
    en = st +1;
    init_Ck = initavg(st:en,:);
    init_Ck_minus = initavg;
    init_Ck_minus(st:en,:) = [];
end
% init_Ck_minus
init_avg_err = init_Ck_minus- kron(ones(m-1,1), init_Ck);
% init_avg_err


%% 
[Kext_k1, Kext_minus_k1]  = gain_gen (Kext, Ck);
 [A_k_k, Q_k] = external_cost_matrices(nCk, nx, Ck, a, b2, Lext, init, Kext_k1, Kext_minus_k1, Rk);

[t_avg, Y_k] = ode23s(@(t,y) A_k_k*y, tspan, init_avg_err);

%%  Separate the components of Y_k into Y_k1 and Y_k2
l = size(Y_k,2);
Y_k1 = [];  % first compnent of x_int
Y_k2 = [];    % second compnent of x_int

for i=  2:2:l
    Y_k1 =[Y_k1 Y_k(:,i-1)];
    Y_k2 =[Y_k2 Y_k(:,i)];
end


%% Instantaneous Cost

avg_cost = zeros(length(Y_k),1);
for i = 1:length(Y_k)
    cost = Y_k(i,:)*Q_k*Y_k(i,:)';
    avg_cost(i) = delta_t*cost;
end

%% Cummulative cost
J_avg_inst  = avg_cost;
% J_cost_avg_cum = cumsum(J_cost_avg);
end












